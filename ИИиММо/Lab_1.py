import os
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from pathlib import Path

# Показывать все столбцы при выводе DataFrame
pd.set_option('display.max_columns', None)

# Путь к корню проекта и датасету
ROOT_DIR = Path(__file__).resolve().parents[1]
file_path = ROOT_DIR / "train.csv"

# --- Проверка существования и корректности файла ---
if not file_path.exists():
    print("Ошибка: файл train.csv не найден.")
    exit()

if os.path.getsize(file_path) == 0:
    print("Ошибка: файл train.csv пустой.")
    exit()

try:
    df = pd.read_csv(file_path)
except pd.errors.EmptyDataError:
    print("Ошибка: файл train.csv пустой или повреждён.")
    exit()
except Exception as e:
    print("Ошибка при чтении файла:", e)
    exit()

# --- Первичный осмотр данных ---
print("=== Первые 5 строк ===")
print(df.head())
print()

print("=== Информация о датасете ===")
df.info()
print()

print("=== Пропуски ДО заполнения ===")
print(df.isnull().sum())
print()

# --- Заполнение пропусков ---

# Возраст заполняем медианой — более устойчива к выбросам
if "Age" in df.columns:
    df["Age"] = df["Age"].fillna(df["Age"].median())

# Числовые признаки расходов заполняем средним значением
expense_cols = ["RoomService", "FoodCourt", "ShoppingMall", "Spa", "VRDeck"]
for col in expense_cols:
    if col in df.columns:
        df[col] = df[col].fillna(df[col].mean())

# Категориальные столбцы заполняем модой (самым частым значением)
cat_fill = ["HomePlanet", "CryoSleep", "Cabin", "Destination", "VIP", "Name"]
for col in cat_fill:
    if col in df.columns and not df[col].mode().empty:
        df[col] = df[col].fillna(df[col].mode()[0])

print("=== Пропуски ПОСЛЕ заполнения ===")
print(df.isnull().sum())
print()

# PassengerId и Name — идентификаторы, Cabin — слишком мало структуры
drop_cols = [c for c in ["PassengerId", "Name", "Cabin"] if c in df.columns]
if drop_cols:
    df = df.drop(columns=drop_cols)

print("=== Столбцы после очистки ===")
print(df.columns.tolist())
print()

num_cols = df.select_dtypes(include=["number"]).columns.tolist()

# Целевую переменную Transported не нормализуем
if "Transported" in num_cols:
    num_cols.remove("Transported")

# MinMaxScaler приводит значения к диапазону [0, 1]
if len(num_cols) > 0:
    scaler = MinMaxScaler()
    df[num_cols] = scaler.fit_transform(df[num_cols])

print("=== После нормализации ===")
print(df.head())
print()

# --- One-hot кодирование категориальных столбцов ---
cat_cols = df.select_dtypes(include=["object", "category", "bool"]).columns.tolist()
if len(cat_cols) > 0:
    # drop_first=True убирает одну фиктивную переменную для избежания мультиколлинеарности
    df = pd.get_dummies(df, columns=cat_cols, drop_first=True)

# Приводим Transported к числовому типу (0/1)
if "Transported" in df.columns:
    df["Transported"] = df["Transported"].astype(int)

print("=== После кодирования ===")
print(df.head())
print()

# --- Сохранение обработанных данных ---
output_path = ROOT_DIR / "lab1_processed_train.csv"
df.to_csv(output_path, index=False)

print(f"Готово! Сохранено {output_path.name}, размер: {df.shape}")
