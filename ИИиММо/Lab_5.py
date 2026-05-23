import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, ConfusionMatrixDisplay

ROOT_DIR = Path(__file__).resolve().parents[1]


def generate_data():
    """Генерация синтетического датасета 'Выборы президента'.
    12 бинарных признаков, 2 класса."""
    np.random.seed(42)

    # 500 объектов, 12 бинарных ответов избирателя (0 или 1)
    X = np.random.randint(0, 2, size=(500, 12))

    # Целевая переменная: правящая партия выигрывает,
    # если сумма первых 6 ответов >= суммы последних 6
    y = (X[:, :6].sum(axis=1) >= X[:, 6:].sum(axis=1)).astype(int)

    # Сохраняем в файлы (требование задания)
    np.savetxt(ROOT_DIR / "lab5_dataIn.txt", X, fmt="%d")
    np.savetxt(ROOT_DIR / "lab5_dataOut.txt", y.reshape(-1, 1), fmt="%d")

    # Считываем обратно из файлов для проверки корректности ввода-вывода
    X_loaded = np.loadtxt(ROOT_DIR / "lab5_dataIn.txt")
    y_loaded = np.loadtxt(ROOT_DIR / "lab5_dataOut.txt").astype(int)

    # Если один объект — убедиться, что матрица двумерная
    if X_loaded.ndim == 1:
        X_loaded = X_loaded.reshape(1, -1)

    return X_loaded, y_loaded


def main():
    X, y = generate_data()

    print(f"Размер X: {X.shape}")
    print(f"Размер y: {y.shape}")
    print(f"Классы: {np.bincount(y)}")

    # stratify=y — сохраняем баланс классов при разделении
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y,
    )

    # Масштабирование признаков (важно для нейросетей и лог. регрессии)
    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc = scaler.transform(X_test)

    # ===== MLP — многослойный перцептрон =====
    mlp = MLPClassifier(
        hidden_layer_sizes=(16,),    # один скрытый слой из 16 нейронов
        activation="logistic",       # сигмоидная функция активации
        solver="adam",               # оптимизатор Adam
        max_iter=1000,
        random_state=42,
    )
    mlp.fit(X_train_sc, y_train)
    y_pred_mlp = mlp.predict(X_test_sc)
    mlp_acc = accuracy_score(y_test, y_pred_mlp)

    print(f"\nMLP Accuracy: {mlp_acc:.4f}")
    print("\nОтчёт MLP:")
    print(classification_report(y_test, y_pred_mlp, digits=4))

    # --- Матрица ошибок ---
    cm = confusion_matrix(y_test, y_pred_mlp)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                                  display_labels=["Правящая партия", "Оппозиция"])
    disp.plot()
    plt.title("Матрица ошибок (MLP)")
    plt.tight_layout()
    plt.savefig(ROOT_DIR / "lab5_confusion_matrix.png", dpi=150)
    plt.close()

    # --- График функции потерь (loss curve) по эпохам ---
    plt.figure(figsize=(8, 5))
    plt.plot(mlp.loss_curve_)
    plt.xlabel("Итерация")
    plt.ylabel("Loss")
    plt.title("Обучение MLP: функция потерь")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(ROOT_DIR / "lab5_loss_curve.png", dpi=150)
    plt.close()

    # ===== Сравнение с другими моделями =====

    # Логистическая регрессия (линейная модель)
    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr.fit(X_train_sc, y_train)
    lr_acc = accuracy_score(y_test, lr.predict(X_test_sc))

    # Random Forest (ансамбль деревьев, не требует масштабирования)
    rf = RandomForestClassifier(n_estimators=200, random_state=42)
    rf.fit(X_train, y_train)
    rf_acc = accuracy_score(y_test, rf.predict(X_test))

    print(f"Logistic Regression Accuracy: {lr_acc:.4f}")
    print(f"Random Forest Accuracy:       {rf_acc:.4f}")

    # --- Столбчатая диаграмма сравнения точности моделей ---
    names = ["MLP", "LogReg", "RandomForest"]
    scores = [mlp_acc, lr_acc, rf_acc]

    plt.figure(figsize=(8, 5))
    bars = plt.bar(names, scores, color=["#4c72b0", "#55a868", "#c44e52"])
    plt.ylim(0, 1)
    # Подписи значений над столбцами
    for bar, score in zip(bars, scores):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                 f"{score:.3f}", ha="center", fontsize=11)
    plt.xlabel("Модель")
    plt.ylabel("Accuracy")
    plt.title("Сравнение моделей")
    plt.tight_layout()
    plt.savefig(ROOT_DIR / "lab5_model_comparison.png", dpi=150)
    plt.close()

    print("Графики сохранены.")


if __name__ == "__main__":
    main()
