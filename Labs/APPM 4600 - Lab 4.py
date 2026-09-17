import matplotlib.pyplot as plt


columns = ["Method", "Input", "Iteration", "Idea Behind Method", "Required for Convergence", "Pros", "Cons"]
rows = [["Bisection", "", "", "", "", "", ""], 
        ["Fixed Point", "", "", "", "", "", ""],
        ["Newton", "", "", "", "", "", ""],
        ["Secant", "", "", "", "", "", ""]]

fig, ax = plt.subplots(figsize=(12,6))

ax.set_axis_off()

table = ax.table(
    cellText=rows,
    colLabels=columns,
    cellLoc="center",
    loc="center",
)


# Increase row height
for row in range(5):
    for col in range(7):
        table[(row, col)]

plt.title("Chart for root finding", fontsize=14)

plt.show()