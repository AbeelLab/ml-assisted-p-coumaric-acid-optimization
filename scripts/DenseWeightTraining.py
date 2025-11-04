""" runs the ML modeling using DenseWeight"""

import numpy as np
from denseweight import DenseWeight
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sp
from sympy import *
import xgboost as xgb
from xgboost import XGBClassifier, XGBRegressor
from sklearn.model_selection import train_test_split
from scipy.stats import linregress

def gradient(predt: np.ndarray,
             dtrain: xgb.DMatrix) -> np.ndarray:
    """Custom gradient based on denseweight balancing for XGboost"""
    y = dtrain.get_label()
    grad = 2 * (predt - y) * dw(y)
    return grad


def hessian(predt: np.ndarray,
            dtrain: xgb.DMatrix) -> np.ndarray:
    "Custom hessian based on  denseweight balancing for XGboost"
    y = dtrain.get_label()
    hessian = 2 * dw(y)
    return hessian


def balanced_mse_wrapper(dw):
    """Wrapper for the balanced mse"""

    def balanced_mse(predt: np.ndarray, dtrain: xgb.DMatrix) -> np.ndarray:
        y = dtrain.get_label()

        grad = 2 * (predt - y) * dw(y)
        hess = 2 * dw(y)
        return grad, hess

    return balanced_mse


# screening results 1 are from the first dbtl round, the remeasured top 100 strains are included

screening_results_1 = pd.read_csv(
    "data/processed/IntegratedData_WURTUD/batch_corrected_screening_results_integrated.csv")
screening_results_2 = pd.read_csv(
    "data/processed/IntegratedData_WURTUD/2504_topX_corrected_screening_results.csv")  #remeasured top 100

screening_results = screening_results_1['batch_corrected_Coumaric_acid'].to_list() + screening_results_2[
    'batch_corrected_Coumaric_acid'].to_list()

# Fit denseweight object
dw = DenseWeight(alpha=0.5, eps=0.01)
density = dw.fit(screening_results)

fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(screening_results, density, s=1)
ax.set_xlabel("[p-CA]")
ax.set_ylabel("Loss weights")
ax.set_title("DenseWeight weighing function")
ax.set_ylim(0, 2.25)
fig.savefig("figures/MachineLearning/DenseWeight_alpha1_eps0.1.png", bbox_inches='tight')
fig.savefig("figures/MachineLearning/DenseWeight_alpha1_eps0.1.svg", bbox_inches='tight')

fig, ax = plt.subplots(figsize=(4, 4))
pd.Series(screening_results).plot(kind="kde")
ax.set_xlabel("[p-CA]")
ax.set_ylabel("Density")

fig.savefig("figures/MachineLearning/DenseWeight_histogram_alpha1_eps0.1.png", bbox_inches='tight')
fig.savefig("figures/MachineLearning/DenseWeight_histogram_alpha1_eps0.1.svg", bbox_inches='tight')

# Denseweight
ypred, ytarget, alpha = sp.symbols('y_{pred} y_{target} alpha')
f = Function('f')(alpha, ytarget)
expr = f * (ypred - ytarget) ** 2
gradient = expr.diff(ypred).simplify()
hessian = gradient.diff(ypred).simplify()

# denseweight
dw = DenseWeight(alpha=0.3, eps=0.1)
weights = dw.fit(screening_results)

# Load data and processing and merging

# 2203_integrated_corrected_screening_results""
# 2504_topX_corrected_screening_results
filename = "data/processed/IntegratedData_WURTUD/strain_numeric_matrix_TEST.csv"
filename_screening1 = ("data/processed/IntegratedData_WURTUD/"
                       "2203_integrated_corrected_screening_TEST.csv")
filename_screening2 = ("data/processed/IntegratedData_WURTUD/"
                       "2504_topX_corrected_screening_results.csv")

strain_count_matrix = pd.read_csv(filename, index_col=0)

strain_count_matrix = strain_count_matrix.iloc[np.where(strain_count_matrix.sum(axis=1) > 0)[0], :]

#screening first round
screening1 = pd.read_csv(filename_screening1, index_col=0)
screening1 = screening1['batch_corrected_Coumaric_acid']

#screening remeasured round
screening2 = pd.read_csv(filename_screening2, index_col=1)
screening2 = screening2['batch_corrected_Coumaric_acid']

screening = pd.concat([screening1, screening2])

data = pd.merge(strain_count_matrix, screening, how="left", left_index=True, right_on="SampleName")

#%%
#ML modelling
train_x, test_x, train_y, test_y = train_test_split(data.iloc[:, :-1], data.iloc[:, -1], test_size=0.20, shuffle=True)

predicted_values = []
true_values = []
r2_values = []
predicted_values_top = []
true_values_top = []
r2_values_top = []

parameters = {'tree_method': 'auto', 'reg_lambda': 1, 'max_depth': 2, "disable_default_eval_metric": 0}
# perform 100 train_test splits, with 10% test set size
for i in range(100):
    train_x, test_x, train_y, test_y = train_test_split(data.iloc[:, :-1], data.iloc[:, -1], test_size=0.20,
                                                        shuffle=True)
    top_test_y_indices = np.argsort(test_y).values[::-1][0:20]

    test_y_top = test_y.iloc[list(top_test_y_indices)]
    test_x_top = test_x.iloc[list(top_test_y_indices), :]

    dtrain = xgb.DMatrix(train_x, label=train_y)
    dtest = xgb.DMatrix(test_x, label=test_y.to_list())

    dtest_top = xgb.DMatrix(test_x_top, label=test_y_top.to_list())

    # a test of the
    regressor = TabPFNRegressor()

    bst = xgb.train(params=parameters,  # any other tree method is fine.
                    dtrain=dtrain,
                    evals=[(dtrain, "train"), (dtest, "validation")],
                    num_boost_round=300,
                    early_stopping_rounds=40,
                    obj=balanced_mse_wrapper(dw))

    preds = bst.predict(dtest)
    true_value = dtest.get_label()

    preds_top = bst.predict(dtest_top)
    true_value_top = dtest_top.get_label()

    slope, intercept, r_value, p_value, std_err = linregress(true_value, preds)
    slope, intercept, r_value_top, p_value, std_err = linregress(true_value_top, preds_top)

    score_GradBoost = r_value ** 2
    score_GradBoost_top = r_value_top ** 2

    r2_values.append(score_GradBoost)
    r2_values_top.append(score_GradBoost_top)
    # print(r2_values)

    predicted_values.append(preds)
    true_values.append(true_value)

    predicted_values_top.append(preds_top)
    true_values_top.append(true_value_top)

#%%


fig, ax = plt.subplots(figsize=(4, 4))
cmap = plt.colormaps["tab20"]
for i in range(5):
    ax.scatter(true_values[i], predicted_values[i], color="#1B1919FF",s=8)

ax.set_xlabel("True [p-CA]")
ax.set_ylabel("Predicted [p-CA]")
ax.set_xlim(0.0, np.max(data.iloc[:, -1] + 0.1))
ax.set_ylim(0.0, np.max(data.iloc[:, -1] + 0.1))
annotate_string = "$R^2$: " + str(np.round(np.mean(r2_values), decimals=2))
annotate_string2 = "R^2 top20: " + str(np.round(np.mean(r2_values_top), decimals=2))

ax.plot([0, 3.3], [0, 3.3], transform=ax.transAxes,linestyle="--",c="black")

ax.text(.04, .96, annotate_string,
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=12)
# ax.text(.01,.94,annotate_string2,ha="left",va="top",transform=ax.transAxes)
plt.title("Predicted versus Measured [p-CA]")
plt.tight_layout()
plt.show()
fig.savefig("figures/MachineLearning/PearsonCorrelation.png",bbox_inches='tight')
fig.savefig("figures/MachineLearning/PearsonCorrelation.svg",bbox_inches='tight')


#%%
# bst.save_model("140624_XGBoost_retrained.model")