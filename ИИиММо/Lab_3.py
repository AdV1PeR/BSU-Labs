from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier, plot_tree
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import accuracy_score, roc_curve, auc
import matplotlib.pyplot as plt

# Загрузка предобработанного датасета из Lab 1
ROOT_DIR = Path(__file__).resolve().parents[1]
DATASET_PATH = ROOT_DIR / "lab1_processed_train.csv"
df = pd.read_csv(DATASET_PATH)


# ============================================================
# 1. РЕГРЕССИЯ ДЕРЕВОМ РЕШЕНИЙ — предсказание Age
# ============================================================

y_reg = df["Age"]
X_reg = df.drop(columns=["Age"])

# Разделение 80/20
X_reg_train, X_reg_test, y_reg_train, y_reg_test = train_test_split(
    X_reg, y_reg, test_size=0.2, random_state=42,
)

# max_depth=6 ограничивает глубину дерева, чтобы снизить переобучение
reg_tree = DecisionTreeRegressor(max_depth=6, random_state=42)
reg_tree.fit(X_reg_train, y_reg_train)
y_reg_pred = reg_tree.predict(X_reg_test)

# Метрики регрессии
reg_mse = mean_squared_error(y_reg_test, y_reg_pred)
reg_rmse = np.sqrt(reg_mse)
reg_mae = mean_absolute_error(y_reg_test, y_reg_pred)
reg_r2 = reg_tree.score(X_reg_test, y_reg_test)  # Коэффициент детерминации


# ============================================================
# 2. КЛАССИФИКАЦИЯ ДЕРЕВОМ РЕШЕНИЙ — предсказание Transported
# ============================================================

target_clf = "Transported_True"
y_clf = df[target_clf].astype(int)
X_clf = df.drop(columns=[target_clf])

X_clf_train, X_clf_test, y_clf_train, y_clf_test = train_test_split(
    X_clf, y_clf, test_size=0.2, random_state=42,
)

clf_tree = DecisionTreeClassifier(max_depth=6, random_state=42)
clf_tree.fit(X_clf_train, y_clf_train)
y_clf_pred = clf_tree.predict(X_clf_test)

# Вероятности положительного класса — нужны для построения ROC-кривой
y_clf_proba = clf_tree.predict_proba(X_clf_test)
y_positive_proba = y_clf_proba[:, 1]

clf_acc = accuracy_score(y_clf_test, y_clf_pred)

# Построение ROC-кривой: зависимость TPR от FPR при разных порогах
fpr, tpr, _ = roc_curve(y_clf_test, y_positive_proba)
roc_auc = auc(fpr, tpr)


# ============================================================
# 3. ГРАФИКИ
# ============================================================

# --- ROC-кривая ---
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, marker="o", markersize=3, label=f"ROC-AUC = {roc_auc:.4f}")
plt.plot([0, 1], [0, 1], linestyle="--", color="grey", label="Случайный классификатор")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC-кривая (дерево решений)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(ROOT_DIR / "lab3_roc_curve.png", dpi=150)
plt.close()

# --- Визуализация дерева регрессии (показываем верхние 3 уровня) ---
plt.figure(figsize=(20, 10))
plot_tree(reg_tree, feature_names=X_reg.columns, filled=True, max_depth=3, fontsize=8)
plt.title("Дерево регрессии (Age)")
plt.tight_layout()
plt.savefig(ROOT_DIR / "lab3_regression_tree.png", dpi=150)
plt.close()

# --- Визуализация дерева классификации ---
plt.figure(figsize=(20, 10))
plot_tree(clf_tree, feature_names=X_clf.columns, class_names=["Not", "Transported"],
          filled=True, max_depth=3, fontsize=8)
plt.title("Дерево классификации (Transported)")
plt.tight_layout()
plt.savefig(ROOT_DIR / "lab3_classification_tree.png", dpi=150)
plt.close()


# --- ВЫВОД РЕЗУЛЬТАТОВ ---
print("=== РЕГРЕССИЯ (Decision Tree) ===")
print(f"MSE:  {reg_mse:.6f}")
print(f"RMSE: {reg_rmse:.6f}")
print(f"MAE:  {reg_mae:.6f}")
print(f"R²:   {reg_r2:.6f}")
print()
print("=== КЛАССИФИКАЦИЯ (Decision Tree) ===")
print(f"Accuracy: {clf_acc:.6f}")
print(f"ROC-AUC:  {roc_auc:.6f}")
