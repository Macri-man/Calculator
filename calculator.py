import streamlit as st
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

st.set_page_config(page_title="Multi-Function Calculator", layout="wide")

st.title("🧮 Multi-Function Calculator")

# Sidebar operation selection
operation = st.sidebar.selectbox("Choose Operation", [
    "Basic Arithmetic",
    "Exponentiation & Logs",
    "Calculus",
    "Statistics",
    "Matrix Operations",
    "Complex Arithmetic",
    "Plot Function",
    "Prime Checker/Generator",
    "Solve System of Equations"
])

# Clear instructions sidebar
with st.sidebar.expander("Instructions"):
    st.write("""
    - Select an operation from above.
    - Enter inputs in the main panel.
    - Results will be shown below inputs.
    - Use valid math expressions (e.g. x**2, sin(x), 3+4j).
    """)

# Helper function for safe sympy eval
def safe_eval(expr):
    try:
        return sp.sympify(expr)
    except:
        st.error(f"Invalid expression: {expr}")
        return None

# --------- Basic Arithmetic -----------
if operation == "Basic Arithmetic":
    st.header("Basic Arithmetic")
    col1, col2 = st.columns(2)
    with col1:
        expr1 = st.text_input("Expression 1", value="10")
    with col2:
        expr2 = st.text_input("Expression 2", value="5")

    op = st.radio("Operation", ["Add", "Subtract", "Multiply", "Divide"])

    val1 = safe_eval(expr1)
    val2 = safe_eval(expr2)
    
    if val1 is not None and val2 is not None:
        result = None
        if op == "Add":
            result = val1 + val2
        elif op == "Subtract":
            result = val1 - val2
        elif op == "Multiply":
            result = val1 * val2
        elif op == "Divide":
            if val2 == 0:
                st.error("Cannot divide by zero!")
            else:
                result = val1 / val2
        if result is not None:
            st.success(f"Result: {result}")

# --------- Exponentiation & Logs ----------
elif operation == "Exponentiation & Logs":
    st.header("Exponentiation & Logarithms")
    sub_op = st.radio("Select Operation", ["Exponentiation", "Logarithm", "Square Root"])

    if sub_op == "Exponentiation":
        base = st.text_input("Base", value="2")
        exponent = st.text_input("Exponent", value="3")
        b = safe_eval(base)
        e = safe_eval(exponent)
        if b is not None and e is not None:
            st.success(f"Result: {b**e}")

    elif sub_op == "Logarithm":
        val = st.text_input("Value", value="10")
        base = st.text_input("Base (default e)", value="e")
        try:
            v = float(val)
            if v <= 0:
                st.error("Value must be positive")
            else:
                if base.lower() == 'e':
                    res = sp.log(v)
                else:
                    b = float(base)
                    if b <= 0 or b == 1:
                        st.error("Base must be positive and not 1")
                    else:
                        res = sp.log(v, b)
                st.success(f"Result: {res.evalf()}")
        except Exception as e:
            st.error(f"Error: {e}")

    else:  # Square Root
        val = st.text_input("Value", value="9")
        try:
            v = float(val)
            if v < 0:
                st.error("Value cannot be negative")
            else:
                st.success(f"Result: {np.sqrt(v)}")
        except Exception as e:
            st.error(f"Error: {e}")

# --------- Calculus ----------
elif operation == "Calculus":
    st.header("Calculus Operations")
    calc_op = st.radio("Operation", ["Derivative", "Nth Derivative", "Integral (Indefinite)", "Integral (Definite)"])
    expr = st.text_input("Function in x", value="x**3 + 2*x")

    x = sp.symbols('x')
    f = safe_eval(expr)

    if f is not None:
        if calc_op == "Derivative":
            d = sp.diff(f, x)
            st.success(f"Derivative: {d}")
        elif calc_op == "Nth Derivative":
            n = st.number_input("Derivative order (n)", min_value=1, value=1)
            d = sp.diff(f, x, int(n))
            st.success(f"{n}th Derivative: {d}")
        elif calc_op == "Integral (Indefinite)":
            integ = sp.integrate(f, x)
            st.success(f"Indefinite Integral: {integ} + C")
        else:
            a = st.text_input("Lower limit (a)", value="0")
            b = st.text_input("Upper limit (b)", value="1")
            try:
                a_val = float(sp.sympify(a))
                b_val = float(sp.sympify(b))
                integ = sp.integrate(f, (x, a_val, b_val))
                st.success(f"Definite Integral from {a_val} to {b_val}: {integ.evalf()}")
            except Exception as e:
                st.error(f"Invalid limits: {e}")

# --------- Statistics ----------
elif operation == "Statistics":
    st.header("Statistics")
    data_str = st.text_area("Enter numbers separated by commas", value="1, 2, 3, 4, 5")
    try:
        data = [float(n.strip()) for n in data_str.split(",") if n.strip()]
        arr = np.array(data)
        st.write(f"Mean: {np.mean(arr)}")
        st.write(f"Median: {np.median(arr)}")
        st.write(f"Variance: {np.var(arr, ddof=1) if len(arr) > 1 else 0}")
        st.write(f"Standard Deviation: {np.std(arr, ddof=1) if len(arr) > 1 else 0}")
        st.write(f"Range: {np.ptp(arr)}")
        counts = {val:data.count(val) for val in set(data)}
        mode = max(counts, key=counts.get) if counts else None
        st.write(f"Mode: {mode}")
    except Exception as e:
        st.error(f"Error parsing numbers: {e}")

# --------- Matrix Operations ----------
elif operation == "Matrix Operations":
    st.header("Matrix Operations")
    mat1_str = st.text_area("Matrix 1 (rows separated by newlines, values by spaces)", value="1 2\n3 4")
    mat2_str = st.text_area("Matrix 2 (optional for Add/Multiply)", value="5 6\n7 8")
    mat_op = st.selectbox("Select operation", ["Add", "Multiply", "Transpose (Matrix 1)", "Determinant (Matrix 1)", "Inverse (Matrix 1)"])

    def parse_matrix(text):
        return np.array([[float(x) for x in row.split()] for row in text.strip().split("\n")])

    try:
        mat1 = parse_matrix(mat1_str)
        mat2 = parse_matrix(mat2_str) if mat2_str.strip() else None

        if mat_op == "Add":
            if mat2 is None:
                st.error("Matrix 2 needed for addition")
            else:
                st.write(mat1 + mat2)
        elif mat_op == "Multiply":
            if mat2 is None:
                st.error("Matrix 2 needed for multiplication")
            else:
                st.write(np.dot(mat1, mat2))
        elif mat_op == "Transpose (Matrix 1)":
            st.write(mat1.T)
        elif mat_op == "Determinant (Matrix 1)":
            st.write(np.linalg.det(mat1))
        elif mat_op == "Inverse (Matrix 1)":
            det = np.linalg.det(mat1)
            if det == 0:
                st.error("Matrix is singular, no inverse")
            else:
                st.write(np.linalg.inv(mat1))
    except Exception as e:
        st.error(f"Matrix input error: {e}")

# --------- Complex Arithmetic ----------
elif operation == "Complex Arithmetic":
    st.header("Complex Number Arithmetic")
    c1 = st.text_input("Complex number 1", value="1+2j")
    c2 = st.text_input("Complex number 2", value="3-4j")
    comp_op = st.selectbox("Operation", ["Add", "Subtract", "Multiply", "Divide", "Exponentiate"])

    try:
        num1 = complex(c1)
        num2 = complex(c2)
        result = None
        if comp_op == "Add":
            result = num1 + num2
        elif comp_op == "Subtract":
            result = num1 - num2
        elif comp_op == "Multiply":
            result = num1 * num2
        elif comp_op == "Divide":
            if num2 == 0:
                st.error("Division by zero")
            else:
                result = num1 / num2
        elif comp_op == "Exponentiate":
            result = num1 ** num2
        if result is not None:
            st.success(f"Result: {result}")
    except Exception as e:
        st.error(f"Error: {e}")

# --------- Plot Function ----------
elif operation == "Plot Function":
    st.header("Function Plotter")
    expr = st.text_input("Enter function of x:", value="sin(x)")
    x_min = st.number_input("x min:", value=-10.0)
    x_max = st.number_input("x max:", value=10.0)
    points = st.slider("Number of points", min_value=100, max_value=1000, value=500)

    try:
        x = sp.symbols('x')
        f = sp.sympify(expr)
        f_np = sp.lambdify(x, f, "numpy")

        x_vals = np.linspace(x_min, x_max, points)
        y_vals = f_np(x_vals)

        fig, ax = plt.subplots()
        ax.plot(x_vals, y_vals)
        ax.set_title(f"Plot of {expr}")
        ax.grid(True)
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Error plotting function: {e}")

# --------- Prime Checker/Generator ----------
elif operation == "Prime Checker/Generator":
    st.header("Prime Checker & Generator")
    choice = st.radio("Select mode", ["Check Prime", "Generate Primes up to N"])

    def is_prime(n):
        if n <= 1:
            return False
        if n <= 3:
            return True
        if n % 2 == 0 or n % 3 == 0:
            return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True

    if choice == "Check Prime":
        num = st.number_input("Enter integer", min_value=1, step=1)
        if st.button("Check"):
            if is_prime(num):
                st.success(f"{num} is prime!")
            else:
                st.warning(f"{num} is not prime.")
    else:
        limit = st.number_input("Generate primes up to", min_value=2, step=1)
        if st.button("Generate"):
            primes = [i for i in range(2, limit + 1) if is_prime(i)]
            st.write(f"Primes up to {limit}:")
            st.write(primes)

# --------- Solve System of Equations ----------
elif operation == "Solve System of Equations":
    st.header("Solve System of Linear Equations")
    eqns_text = st.text_area("Enter equations (one per line, e.g. x + y = 2)", value="x + y = 2\nx - y = 0")
    vars_text = st.text_input("Variables (comma separated)", value="x,y")

    try:
        eqns = [sp.Eq(*map(sp.sympify, eq.replace("=", "-(").split("-("))) for eq in eqns_text.strip().split("\n")]
        variables = sp.symbols(vars_text)
        sol = sp.solve(eqns, variables)
        st.write("Solution:")
        st.write(sol)
    except Exception as e:
        st.error(f"Error solving system: {e}")
