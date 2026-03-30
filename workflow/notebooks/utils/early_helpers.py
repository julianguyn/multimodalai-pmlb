import pandas as pd
from sklearn.model_selection import KFold, train_test_split, GridSearchCV
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor, AdaBoostRegressor
from sklearn.metrics import root_mean_squared_error, mean_absolute_error, make_scorer
import xgboost as xgb
from scipy.stats import spearmanr, pearsonr
import shap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
newcmp = mcolors.LinearSegmentedColormap.from_list("custom", ["blue", "red"])


def spearman_scorer(y_true, y_pred):
    """
    Helper function to create scorer
    """
    return spearmanr(y_true, y_pred).correlation


def run_shap(X, final_model, filepath):
    """
    Evaluate SHAP on selected model
    """

    # run shap
    explainer = shap.Explainer(final_model, X)
    shap_values = explainer(X)

    # summary plot
    plt.figure()
    shap.summary_plot(shap_values, X, show=False)
    for fc in plt.gcf().get_children():
        for fcc in fc.get_children():
            if hasattr(fcc, "set_cmap"):
                fcc.set_cmap(newcmp)
    plt.savefig(("data/results/figures/"+filepath+"shap_summary.png"), dpi=300, bbox_inches='tight')
    plt.close()

    # heatmap
    plt.figure()
    shap.plots.heatmap(
        shap_values,
        instance_order=shap_values.sum(1),
        show=False
    )
    plt.savefig(("data/results/figures/"+filepath+"heatmap.png"), dpi=300, bbox_inches='tight')
    plt.close()


"""
Simple ML models with nested 5foldCV
    * LASSO
    * Elastic Net
    * Random Forest
"""

def run_lasso(X, y, path):

    # initialize outer folds (5 folds, 80% train, 20% test)
    outer_cv = KFold(n_splits=5, shuffle=True, random_state=101)

    # initialize variables to store results
    results = []
    hyperparams_list = []
    fold = 1

    # loop through each of the outer five folds
    for train_index, test_index in outer_cv.split(X):

        # initialize inner folds (5 folds, 80% train, 20% test)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=fold)

        # split train and test
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        # initialize LASSO model
        lasso = linear_model.Lasso()

        # specify parameters for optimization
        parameters = {
            'alpha': [0.001, 0.01, 0.1, 1, 10, 100],
            'max_iter': [500, 1000, 5000, 7500]
        }

        # identify optimal parameters
        reg = GridSearchCV(
            estimator = lasso,
            param_grid = parameters,
            cv=inner_cv,
            scoring=make_scorer(spearman_scorer, greater_is_better=True),
            n_jobs=-1
        )

        # fit model
        reg.fit(X_train, y_train)

        # get best model & parameters
        reg_best = reg.best_estimator_

        best_params = reg.best_params_
        hyperparams_list.append(best_params)

        # get predicted values for test data
        y_pred = reg_best.predict(X_test)

        # compute metrics
        s_corr = spearmanr(y_test, y_pred).correlation
        p_corr = pearsonr(y_test, y_pred)[0]
        rmse = root_mean_squared_error(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)

        # save model correlation and features
        results.append({
            "fold": fold,
            "spearman": s_corr,
            "pearson": p_corr,
            "rmse": rmse,
            "mae": mae,
            **best_params
        })

        fold += 1

    # save results to dataframe
    results_df = pd.DataFrame(results)
    results_df.to_csv(("data/results/data/"+path+"lasso_folds.csv"), index=False)

    # get most common hyperparameters
    #params_df = pd.DataFrame(hyperparams_list)
    #final_params = params_df.mode().iloc[0].to_dict()
    #final_params["max_iter"] = int(final_params["max_iter"])

    # fit final model
    #final_model = linear_model.Lasso(**final_params)
    #final_model.fit(X, y)

    # run shap
    #run_shap(X, final_model, path)


def run_elastic_net(X, y, path):

    # initialize outer folds (5 folds, 80% train, 20% test)
    outer_cv = KFold(n_splits=5, shuffle=True, random_state=101)

    # initialize variables to store results
    results = []
    hyperparams_list = []
    fold = 1

    # loop through each of the outer five folds
    for train_index, test_index in outer_cv.split(X):

        # initialize inner folds (5 folds, 80% train, 20% test)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=fold)

        # split train and test
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        # initialize ElasticNet model
        en = linear_model.ElasticNet()

        # specify parameters for optimization
        parameters = {
            'alpha': [0.1, 1, 10, 100],
            'l1_ratio': [0.2, 0.5, 0.8],
            'max_iter': [1000, 5000, 7500]
        }

        # identify optimal parameters
        reg = GridSearchCV(
            estimator = en,
            param_grid = parameters,
            cv=inner_cv,
            scoring=make_scorer(spearman_scorer, greater_is_better=True),
            n_jobs=-1
        )

        # fit model
        reg.fit(X_train, y_train)

        # get best model & parameters
        reg_best = reg.best_estimator_

        best_params = reg.best_params_
        hyperparams_list.append(best_params)

        # get predicted values for test data
        y_pred = reg_best.predict(X_test)

        # compute metrics
        s_corr = spearmanr(y_test, y_pred).correlation
        p_corr = pearsonr(y_test, y_pred)[0]
        rmse = root_mean_squared_error(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)

        # save model correlation and features
        results.append({
            "fold": fold,
            "spearman": s_corr,
            "pearson": p_corr,
            "rmse": rmse,
            "mae": mae,
            **best_params
        })

        fold += 1

    # save results to dataframe
    results_df = pd.DataFrame(results)
    results_df.to_csv(("data/results/data/"+path+"en_folds.csv"), index=False)

    # get most common hyperparameters
    #params_df = pd.DataFrame(hyperparams_list)
    #final_params = params_df.mode().iloc[0].to_dict()

    # fit final model
    #final_model = linear_model.ElasticNet(**final_params)
    #final_model.fit(X, y)

    # run shap
    #run_shap(X, final_model, path)

def run_random_forest(X, y, path):

    # initialize outer folds (5 folds, 80% train, 20% test)
    outer_cv = KFold(n_splits=5, shuffle=True, random_state=101)

    # initialize variables to store results
    results = []
    hyperparams_list = []
    fold = 1

    # loop through each of the outer five folds
    for train_index, test_index in outer_cv.split(X):

        # initialize inner folds (5 folds, 80% train, 20% test)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=fold)

        # split train and test
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        # initialize random forest model
        rf = RandomForestRegressor()

        # specify parameters for optimization
        parameters = {
            'n_estimators': [10, 50, 100, 150, 200],
            'max_depth': [None, 10, 20],
            'min_samples_split': [2, 5],
            'min_samples_leaf': [1, 2, 5],
            'max_features': ['sqrt', 'log2']
        }

        # identify optimal parameters
        reg = GridSearchCV(
            estimator = rf,
            param_grid = parameters,
            cv=inner_cv,
            scoring=make_scorer(spearman_scorer, greater_is_better=True),
            n_jobs=-1
        )

        # fit model
        reg.fit(X_train, y_train)

        # get best model & parameters
        reg_best = reg.best_estimator_

        best_params = reg.best_params_
        hyperparams_list.append(best_params)

        # get predicted values for test data
        y_pred = reg_best.predict(X_test)

        # compute metrics
        s_corr = spearmanr(y_test, y_pred).correlation
        p_corr = pearsonr(y_test, y_pred)[0]
        rmse = root_mean_squared_error(y_test, y_pred)
        mae = mean_absolute_error(y_test, y_pred)

        # save model correlation and features
        results.append({
            "fold": fold,
            "spearman": s_corr,
            "pearson": p_corr,
            "rmse": rmse,
            "mae": mae,
            **best_params
        })

        fold += 1

    # save results to dataframe
    results_df = pd.DataFrame(results)
    results_df.to_csv(("data/results/data/"+path+"rf_folds.csv"), index=False)

    # get most common hyperparameters
    #params_df = pd.DataFrame(hyperparams_list)
    #final_params = params_df.mode().iloc[0].to_dict()
    #final_params["max_depth"] = int(final_params["max_depth"])
    #final_params["min_samples_leaf"] = int(final_params["min_samples_leaf"])

    # fit final model
    #final_model = RandomForestRegressor(**final_params)
    #final_model.fit(X, y)

    # run shap
    #run_shap(X, final_model, path)