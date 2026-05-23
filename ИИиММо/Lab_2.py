from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.metrics import precision_score, recall_score, f1_score
import matplotlib.pyplot as plt


# --- Загрузка исходных данных ---
ROOT_DIR = Path(__file__).resolve().parents[1]
DATASET_PATH = ROOT_DIR / "train.csv"
df = pd.read_csv(DATASET_PATH)

# Извлекаем ID группы из PassengerId (формат "XXXX_YY") и считаем размер группы
df["GroupId"] = df["PassengerId"].str.split("_").str[0]
df["GroupSize"] = df.groupby("GroupId")["PassengerId"].transform("count")

# Разбиваем Cabin (формат "Deck/Num/Side") на три отдельных признака
cabin_split = df["Cabin"].str.split("/", expand=True)
df["CabinDeck"] = cabin_split[0]
df["CabinNum"] = pd.to_numeric(cabin_split[1], errors="coerce")
df["CabinSide"] = cabin_split[2]

# Убираем столбцы, не несущие прямой пользы для модели
df = df.drop(columns=["PassengerId", "Cabin", "Name"])


# ============================================================
# 1) ЗАДАЧА РЕГРЕССИИ — предсказание возраста пассажира (Age)
# ============================================================

X_reg = df.drop(columns=["Age"])
y_reg = df["Age"]

# Работаем только со строками, где Age известен
mask_reg = y_reg.notna()
X_reg = X_reg.loc[mask_reg]
y_reg = y_reg.loc[mask_reg]

# Автоматически определяем типы столбцов для разных пайплайнов
cat_cols_reg = X_reg.select_dtypes(include=["object", "bool"]).columns.tolist()
num_cols_reg = X_reg.select_dtypes(include=["number"]).columns.tolist()

# Числовые: заполнение медианой + стандартизация
num_pipe_reg = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
])

# Категориальные: заполнение модой + one-hot кодирование
cat_pipe_reg = Pipeline([
    ("imputer", SimpleImputer(strategy="most_frequent")),
    ("onehot", OneHotEncoder(handle_unknown="ignore")),
])

# Объединяем оба пайплайна в один трансформер
preprocessor_reg = ColumnTransformer([
    ("num", num_pipe_reg, num_cols_reg),
    ("cat", cat_pipe_reg, cat_cols_reg),
])

# Разделение на обучающую и тестовую выборки (80/20)
X_reg_train, X_reg_test, y_reg_train, y_reg_test = train_test_split(
    X_reg, y_reg, test_size=0.2, random_state=42,
)

# --- Линейная регрессия (без регуляризации) ---
lr_pipeline = Pipeline([
    ("prep", preprocessor_reg),
    ("model", LinearRegression()),
])
lr_pipeline.fit(X_reg_train, y_reg_train)
y_lr_pred = lr_pipeline.predict(X_reg_test)

# Метрики качества регрессии
lr_mse = mean_squared_error(y_reg_test, y_lr_pred)
lr_rmse = np.sqrt(lr_mse)
lr_mae = mean_absolute_error(y_reg_test, y_lr_pred)

# --- Ridge-регрессия (L2-регуляризация для борьбы с переобучением) ---
ridge_pipeline = Pipeline([
    ("prep", preprocessor_reg),
    ("model", Ridge(alpha=1.0)),
])
ridge_pipeline.fit(X_reg_train, y_reg_train)
y_ridge_pred = ridge_pipeline.predict(X_reg_test)

ridge_mse = mean_squared_error(y_reg_test, y_ridge_pred)
ridge_rmse = np.sqrt(ridge_mse)
ridge_mae = mean_absolute_error(y_reg_test, y_ridge_pred)


# ============================================================
# 2) ЗАДАЧА КЛАССИФИКАЦИИ — предсказание Transported
# ============================================================

X_clf = df.drop(columns=["Transported"])
y_clf = df["Transported"].astype(int)

cat_cols_clf = X_clf.select_dtypes(include=["object", "bool"]).columns.tolist()
num_cols_clf = X_clf.select_dtypes(include=["number"]).columns.tolist()

# Аналогичные пайплайны для задачи классификации
num_pipe_clf = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
])

cat_pipe_clf = Pipeline([
    ("imputer", SimpleImputer(strategy="most_frequent")),
    ("onehot", OneHotEncoder(handle_unknown="ignore")),
])

preprocessor_clf = ColumnTransformer([
    ("num", num_pipe_clf, num_cols_clf),
    ("cat", cat_pipe_clf, cat_cols_clf),
])

X_clf_train, X_clf_test, y_clf_train, y_clf_test = train_test_split(
    X_clf, y_clf, test_size=0.2, random_state=42,
)

# Логистическая регрессия — базовый метод бинарной классификации
log_pipeline = Pipeline([
    ("prep", preprocessor_clf),
    ("model", LogisticRegression(max_iter=2000)),
])
log_pipeline.fit(X_clf_train, y_clf_train)
y_clf_pred = log_pipeline.predict(X_clf_test)

# Вычисление метрик классификации
clf_acc = accuracy_score(y_clf_test, y_clf_pred)
clf_prec = precision_score(y_clf_test, y_clf_pred)
clf_rec = recall_score(y_clf_test, y_clf_pred)
clf_f1 = f1_score(y_clf_test, y_clf_pred)
clf_cm = confusion_matrix(y_clf_test, y_clf_pred)
clf_report = classification_report(y_clf_test, y_clf_pred)

# --- Визуализация матрицы ошибок ---
fig, ax = plt.subplots(figsize=(5, 4))
im = ax.imshow(clf_cm, interpolation="nearest", cmap="coolwarm")
ax.set_title("Матрица ошибок (Lab 2)")
fig.colorbar(im, ax=ax)
# Отображаем числовые значения в каждой ячейке
for i in range(clf_cm.shape[0]):
    for j in range(clf_cm.shape[1]):
        ax.text(j, i, str(clf_cm[i, j]), ha="center", va="center", fontsize=14)
ax.set_ylabel("Истинное значение")
ax.set_xlabel("Предсказание")
fig.tight_layout()
fig.savefig(ROOT_DIR / "lab2_confusion_matrix.png")
plt.close(fig)


# --- ВЫВОД РЕЗУЛЬТАТОВ ---
print("=== РЕГРЕССИЯ: предсказание Age ===")
print(f"Linear Regression MSE:  {lr_mse:.6f}")
print(f"Linear Regression RMSE: {lr_rmse:.6f}")
print(f"Linear Regression MAE:  {lr_mae:.6f}")
print()
print("=== УЛУЧШЕНИЕ: Ridge ===")
print(f"Ridge MSE:  {ridge_mse:.6f}")
print(f"Ridge RMSE: {ridge_rmse:.6f}")
print(f"Ridge MAE:  {ridge_mae:.6f}")
print()
print("=== КЛАССИФИКАЦИЯ: предсказание Transported ===")
print(f"Accuracy:  {clf_acc:.6f}")
print(f"Precision: {clf_prec:.6f}")
print(f"Recall:    {clf_rec:.6f}")
print(f"F1-score:  {clf_f1:.6f}")
print()
print("Confusion matrix:")
print(clf_cm)
print()
print("Classification report:")
print(clf_report)
print("Файл lab2_confusion_matrix.png сохранён.")
