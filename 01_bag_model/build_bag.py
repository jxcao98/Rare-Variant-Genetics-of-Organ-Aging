#!/usr/bin/env python3
"""Train a linear SVR biological-age model and calculate corrected BAG values."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVR


RANDOM_SEED = 66


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organ", required=True, help="Organ-system label used in output names.")
    parser.add_argument("--input", required=True, type=Path, help="Tab-delimited participant table.")
    parser.add_argument("--features", required=True, type=Path, help="One biomarker column name per line.")
    parser.add_argument("--test-ids", required=True, type=Path, help="One held-out participant ID per line.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--id-column", default="eid")
    parser.add_argument("--age-column", default="Age")
    parser.add_argument("--sex-column", default="Sex")
    parser.add_argument("--healthy-column", default="AllDisease_flag")
    parser.add_argument("--healthy-value", default="0")
    parser.add_argument("--max-iter", type=int, default=1_000_000)
    return parser.parse_args()


def read_feature_names(path: Path) -> list[str]:
    features = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not features:
        raise ValueError(f"No feature names were found in {path}")
    if len(features) != len(set(features)):
        raise ValueError("The feature list contains duplicate column names")
    return features


def inverse_normal_transform(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    observed = values.dropna()
    ranks = observed.rank(method="average")
    probabilities = (ranks - 0.5) / len(observed)
    result.loc[observed.index] = norm.ppf(probabilities)
    return result


def coerce_healthy_flag(series: pd.Series, healthy_value: str) -> pd.Series:
    return series.astype(str).str.strip().eq(str(healthy_value))


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(args.input, sep="\t", low_memory=False)
    biomarkers = read_feature_names(args.features)
    model_features = biomarkers + [args.sex_column]
    required = {
        args.id_column,
        args.age_column,
        args.sex_column,
        args.healthy_column,
        *biomarkers,
    }
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"Input table is missing columns: {', '.join(missing)}")

    test_ids = set(
        pd.read_csv(args.test_ids, header=None, dtype=str).iloc[:, 0].str.strip().dropna()
    )
    if not test_ids:
        raise ValueError("The held-out test-ID file is empty")

    data[args.id_column] = data[args.id_column].astype(str)
    complete = data.dropna(subset=model_features + [args.age_column]).copy()
    healthy = complete[coerce_healthy_flag(complete[args.healthy_column], args.healthy_value)].copy()
    train = healthy[~healthy[args.id_column].isin(test_ids)].copy()
    test = healthy[healthy[args.id_column].isin(test_ids)].copy()
    if train.empty or test.empty:
        raise ValueError("Healthy training and held-out test sets must both contain complete observations")

    x_train = train[model_features].to_numpy(dtype=float)
    y_train = train[args.age_column].to_numpy(dtype=float)
    y_sd = float(np.std(y_train, ddof=1))
    if not np.isfinite(y_sd) or y_sd == 0:
        raise ValueError("Chronological age has zero or undefined variance in the training set")

    pipeline = Pipeline(
        steps=[
            ("scale", StandardScaler()),
            (
                "svr",
                LinearSVR(
                    random_state=RANDOM_SEED,
                    max_iter=args.max_iter,
                ),
            ),
        ]
    )
    parameter_grid = {
        "svr__C": [0.1, 1.0, 10.0],
        "svr__epsilon": [0.01 * y_sd, 0.1 * y_sd, 0.5 * y_sd],
        "svr__tol": [1e-4, 1e-3, 1e-2],
    }
    cross_validation = KFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
    search = GridSearchCV(
        pipeline,
        parameter_grid,
        scoring="neg_mean_absolute_error",
        cv=cross_validation,
        n_jobs=-1,
        refit=True,
    )
    search.fit(x_train, y_train)
    model = search.best_estimator_

    test_predictions = model.predict(test[model_features].to_numpy(dtype=float))
    test_age = test[args.age_column].to_numpy(dtype=float)
    test_r = pearsonr(test_age, test_predictions).statistic if len(test) > 1 else np.nan

    train_predictions = model.predict(x_train)
    train_raw_bag = train_predictions - y_train
    bias_model = LinearRegression().fit(train[[args.age_column]], train_raw_bag)

    complete_predictions = model.predict(complete[model_features].to_numpy(dtype=float))
    raw_bag = complete_predictions - complete[args.age_column].to_numpy(dtype=float)
    predicted_bias = bias_model.predict(complete[[args.age_column]])
    adjusted_bag = raw_bag - predicted_bias

    results = data[[args.id_column, args.age_column, args.sex_column]].copy()
    results["predicted_age"] = np.nan
    results["bag_raw"] = np.nan
    results["bag_adjusted"] = np.nan
    results.loc[complete.index, "predicted_age"] = complete_predictions
    results.loc[complete.index, "bag_raw"] = raw_bag
    results.loc[complete.index, "bag_adjusted"] = adjusted_bag
    results["bag_adjusted_int"] = inverse_normal_transform(results["bag_adjusted"])

    stem = args.organ.lower().replace(" ", "_")
    results.to_csv(args.output_dir / f"{stem}_bag.tsv", sep="\t", index=False)
    joblib.dump(model, args.output_dir / f"{stem}_linear_svr.joblib")
    joblib.dump(bias_model, args.output_dir / f"{stem}_age_bias_model.joblib")

    metadata = {
        "organ": args.organ,
        "random_seed": RANDOM_SEED,
        "features": model_features,
        "healthy_training_n": int(len(train)),
        "healthy_test_n": int(len(test)),
        "prediction_n": int(len(complete)),
        "best_parameters": search.best_params_,
        "cross_validated_mae": float(-search.best_score_),
        "test_mae": float(mean_absolute_error(test_age, test_predictions)),
        "test_pearson_r": float(test_r),
        "bag_definition": "predicted_age - chronological_age",
        "bias_correction": "training-set raw BAG regressed on chronological age",
    }
    (args.output_dir / f"{stem}_model_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
