from typing import List
import numpy as np

def predict_5_steps(A: List[List[float]], x0: List[float], steps: int = 5) -> List[np.ndarray]:
    A = np.array(A, dtype=float)
    x = np.array(x0, dtype=float)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError(f"Matrix A must be square. Got shape {A.shape}.")
    if x.ndim != 1 or x.shape[0] != A.shape[1]:
        raise ValueError(f"Vector x0 length must match A size. Got len(x0)={x.shape[0]}, A is {A.shape}.")

    predictions = []
    for t in range(1, steps + 1):
        x = A @ x
        predictions.append(x.copy())

    return predictions


if __name__ == "__main__":
    A = [
        [0.9, 0.1],
        [0.2, 0.8]
    ]
    x0 = [100, 50]

    preds = predict_5_steps(A, x0, steps=5)

    print("Initial x0:", np.array(x0, dtype=float))
    for i, xt in enumerate(preds, start=1):
        print(f"x{i} =", xt)
