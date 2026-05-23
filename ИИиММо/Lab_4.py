from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import roc_curve, auc, roc_auc_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt

# Загрузка предобработанного датасета
ROOT_DIR = Path(__file__).resolve().parents[1]
DATASET_PATH = ROOT_DIR / "lab1_processed_train.csv"
df = pd.read_csv(DATASET_PATH)

target_col = "Transported_True"

X = df.drop(columns=[target_col]).astype(float)
y = df[target_col].astype(int)

# stratify=y — сохраняем пропорции классов в train/test
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y,
)


# ============================================================
# 1. RANDOM FOREST — бэггинг деревьев решений
# ============================================================

rf = RandomForestClassifier(
    n_estimators=200,   # кол-во деревьев в лесу
    oob_score=True,     # Out-Of-Bag оценка (на данных, не попавших в бутстреп)
    n_jobs=-1,          # использовать все ядра CPU
    random_state=42,
)
rf.fit(X_train, y_train)
rf_pred = rf.predict(X_test)
rf_proba = rf.predict_proba(X_test)[:, 1]       # вероятности для ROC
rf_oob_proba = rf.oob_decision_function_[:, 1]  # OOB-вероятности

# Метрики на тестовой выборке
rf_acc = accuracy_score(y_test, rf_pred)
rf_prec = precision_score(y_test, rf_pred)
rf_rec = recall_score(y_test, rf_pred)
rf_f1 = f1_score(y_test, rf_pred)
rf_auc = roc_auc_score(y_test, rf_proba)

# OOB-метрики — оценка без отдельной валидации
rf_oob_acc = rf.oob_score_
rf_oob_err = 1 - rf_oob_acc
rf_oob_auc = roc_auc_score(y_train, rf_oob_proba)
rf_report = classification_report(y_test, rf_pred)
rf_cm = confusion_matrix(y_test, rf_pred)


# ============================================================
# 2. ADABOOST — последовательный бустинг слабых классификаторов
# ============================================================

# Базовая модель — пень (дерево глубиной 1)
ada_base = DecisionTreeClassifier(max_depth=1, random_state=42)

ada = AdaBoostClassifier(
    estimator=ada_base,
    n_estimators=150,    # кол-во итераций бустинга
    learning_rate=0.5,   # вклад каждого классификатора
    random_state=42,
)
ada.fit(X_train, y_train)
ada_pred = ada.predict(X_test)
ada_proba = ada.predict_proba(X_test)[:, 1]

ada_acc = accuracy_score(y_test, ada_pred)
ada_prec = precision_score(y_test, ada_pred)
ada_rec = recall_score(y_test, ada_pred)
ada_f1 = f1_score(y_test, ada_pred)
ada_auc = roc_auc_score(y_test, ada_proba)
ada_report = classification_report(y_test, ada_pred)
ada_cm = confusion_matrix(y_test, ada_pred)


# ============================================================
# 3. GRADIENT BOOSTING — градиентный бустинг
# ============================================================

gb = GradientBoostingClassifier(
    n_estimators=200,
    learning_rate=0.1,  # шаг обучения (меньше = точнее, но дольше)
    max_depth=3,        # глубина каждого дерева
    random_state=42,
)
gb.fit(X_train, y_train)
gb_pred = gb.predict(X_test)
gb_proba = gb.predict_proba(X_test)[:, 1]

gb_acc = accuracy_score(y_test, gb_pred)
gb_prec = precision_score(y_test, gb_pred)
gb_rec = recall_score(y_test, gb_pred)
gb_f1 = f1_score(y_test, gb_pred)
gb_auc = roc_auc_score(y_test, gb_proba)
gb_report = classification_report(y_test, gb_pred)
gb_cm = confusion_matrix(y_test, gb_pred)


# ============================================================
# 4. ROC-КРИВЫЕ — сравнение всех трёх моделей
# ============================================================

rf_fpr, rf_tpr, _ = roc_curve(y_test, rf_proba)
ada_fpr, ada_tpr, _ = roc_curve(y_test, ada_proba)
gb_fpr, gb_tpr, _ = roc_curve(y_test, gb_proba)

plt.figure(figsize=(8, 6))
plt.plot(rf_fpr, rf_tpr, label=f"Random Forest (AUC = {auc(rf_fpr, rf_tpr):.4f})")
plt.plot(ada_fpr, ada_tpr, label=f"AdaBoost (AUC = {auc(ada_fpr, ada_tpr):.4f})")
plt.plot(gb_fpr, gb_tpr, label=f"Gradient Boosting (AUC = {auc(gb_fpr, gb_tpr):.4f})")
plt.plot([0, 1], [0, 1], linestyle="--", color="grey", label="Случайный классификатор")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC-кривые ансамблевых методов")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(ROOT_DIR / "lab4_roc_curves.png", dpi=150)
plt.close()


# ============================================================
# 5. ВЫВОД РЕЗУЛЬТАТОВ
# ============================================================

print("=== RANDOM FOREST ===")
print(f"Test Accuracy:  {rf_acc:.6f}")
print(f"Test Precision: {rf_prec:.6f}")
print(f"Test Recall:    {rf_rec:.6f}")
print(f"Test F1:        {rf_f1:.6f}")
print(f"Test ROC-AUC:   {rf_auc:.6f}")
print(f"OOB Accuracy:   {rf_oob_acc:.6f}")
print(f"OOB Error:      {rf_oob_err:.6f}")
print(f"OOB ROC-AUC:    {rf_oob_auc:.6f}")
print("Confusion matrix:")
print(rf_cm)
print(rf_report)

print("=== ADABOOST ===")
print(f"Accuracy:  {ada_acc:.6f}")
print(f"Precision: {ada_prec:.6f}")
print(f"Recall:    {ada_rec:.6f}")
print(f"F1:        {ada_f1:.6f}")
print(f"ROC-AUC:   {ada_auc:.6f}")
print("Confusion matrix:")
print(ada_cm)
print(ada_report)

print("=== GRADIENT BOOSTING ===")
print(f"Accuracy:  {gb_acc:.6f}")
print(f"Precision: {gb_prec:.6f}")
print(f"Recall:    {gb_rec:.6f}")
print(f"F1:        {gb_f1:.6f}")
print(f"ROC-AUC:   {gb_auc:.6f}")
print("Confusion matrix:")
print(gb_cm)
print(gb_report)
