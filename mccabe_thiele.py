# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "numpy",
#     "scipy",
#     "matplotlib",
# ]
# ///

import marimo

__generated_with = "0.14.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import brentq

    return brentq, mo, np, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # The McCabe-Thiele Method

    An interactive, step-by-step derivation of binary distillation design.
    Every equation used by the code below is derived here first.

    **How to read this notebook:** Sections 1-7 build the theory. Section 8 implements it
    as short numpy functions. Section 9 is an interactive explorer. Section 10 has exercises.

    Run it with `uvx marimo edit --sandbox mccabe_thiele.py` (dependencies are
    installed from the PEP 723 header automatically).
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 1. Problem statement

    A column separates a binary mixture of a **light** (more volatile) component A
    and a **heavy** component B. All compositions are **mole fractions of A**.

    | Symbol | Meaning |
    |---|---|
    | $F, D, B$ | feed, distillate, bottoms molar flows |
    | $z_F, x_D, x_B$ | feed, distillate, bottoms mole fraction of A |
    | $L, V$ | liquid and vapor flows in the **rectifying** (top) section |
    | $\bar L, \bar V$ | liquid and vapor flows in the **stripping** (bottom) section |
    | $R = L/D$ | reflux ratio |
    | $q$ | feed thermal condition (defined in Section 5) |
    | $\alpha$ | relative volatility |

    **Goal:** given $F, z_F, x_D, x_B, q, R$ and equilibrium data, find the number of
    equilibrium stages and the optimal feed stage.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 2. Vapor-liquid equilibrium

    For an ideal mixture (Raoult's law, ideal vapor), the partial pressure of each
    component equals its vapor pressure times its liquid mole fraction:

    $$
    y_A P = x_A P_A^{sat}, \qquad y_B P = x_B P_B^{sat}
    $$

    The **relative volatility** is the ratio of the K-values $K_i = y_i/x_i$:

    $$
    \alpha = \frac{y_A/x_A}{y_B/x_B} = \frac{P_A^{sat}}{P_B^{sat}}
    $$

    Write $x = x_A$, $y = y_A$ and use $x_B = 1-x$, $y_B = 1-y$:

    $$
    \alpha = \frac{y\,(1-x)}{x\,(1-y)}
    $$

    **Solve for $y$.** Multiply out: $\alpha x - \alpha x y = y - x y$. Collect the terms
    in $y$:

    $$
    \alpha x = y\,\bigl(1 - x + \alpha x\bigr)
    \quad\Longrightarrow\quad
    \boxed{\,y = \frac{\alpha x}{1 + (\alpha - 1)\,x}\,}
    $$

    **Solve for $x$.** Starting from $\alpha x - \alpha x y = y - x y$, collect terms in $x$:

    $$
    x\,\bigl(\alpha - \alpha y + y\bigr) = y
    \quad\Longrightarrow\quad
    \boxed{\,x = \frac{y}{\alpha - (\alpha - 1)\,y}\,}
    $$

    The second form is what we use to step **horizontally** from a vapor composition to the
    liquid in equilibrium with it. For $\alpha > 1$ the curve lies above the diagonal $y = x$.
    Larger $\alpha$ means easier separation.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 3. Constant molal overflow (CMO)

    The method assumes that the molar latent heats of A and B are equal, heat of mixing is
    negligible and the column is adiabatic. Then every mole of vapor condensed
    vaporizes exactly one mole of liquid, so

    $$
    L, V \text{ are constant within a section, and } \bar L, \bar V \text{ are constant in the other.}
    $$

    The flows only change at the feed stage. This turns the energy balances into
    trivial statements, leaving **material balances only**.

    ## 4. Material balances and the operating lines

    ### 4.1 Overall balances

    $$
    F = D + B \qquad\qquad F z_F = D x_D + B x_B
    $$

    Substitute $B = F - D$ into the component balance:

    $$
    F z_F = D x_D + (F - D)\, x_B
    \;\Longrightarrow\;
    \boxed{\frac{D}{F} = \frac{z_F - x_B}{x_D - x_B}}
    \qquad
    \boxed{\frac{B}{F} = \frac{x_D - z_F}{x_D - x_B}}
    $$

    This also shows why a feasible design needs $x_B < z_F < x_D$.

    ### 4.2 Rectifying section (above the feed)

    Draw a control volume around the top of the column: the total condenser plus stages
    $1 \dots n$. Liquid $L$ leaves stage $n$ (composition $x_n$) and vapor $V$ enters
    from stage $n+1$ (composition $y_{n+1}$).

    Total balance: $V = L + D$.

    Component A balance:

    $$
    V\, y_{n+1} = L\, x_n + D\, x_D
    $$

    Divide by $V = L + D$ and use $R = L/D$, so $L/V = R/(R+1)$ and $D/V = 1/(R+1)$:

    $$
    \boxed{\,y_{n+1} = \frac{R}{R+1}\,x_n + \frac{x_D}{R+1}\,}
    $$

    **Checks:** slope is $R/(R+1) < 1$. The intercept at $x=0$ is $x_D/(R+1)$. At $x = x_D$:
    $\frac{R\,x_D + x_D}{R+1} = x_D$, so the line passes through $(x_D, x_D)$ on the diagonal.

    ### 4.3 Stripping section (below the feed)

    Draw a control volume around the bottom: the reboiler plus stages below the feed.
    Liquid $\bar L$ enters from stage $m$ (composition $x_m$) and vapor $\bar V$ leaves the stage
    below it (composition $y_{m+1}$).

    Total balance: $\bar L = \bar V + B$.

    Component A balance:

    $$
    \bar L\, x_m = \bar V\, y_{m+1} + B\, x_B
    $$

    Solve for $y_{m+1}$:

    $$
    \boxed{\,y_{m+1} = \frac{\bar L}{\bar V}\,x_m - \frac{B}{\bar V}\,x_B\,}
    $$

    **Check:** at $x_m = x_B$ we get $y = \frac{(\bar L - B)}{\bar V} x_B = x_B$, so the line passes
    through $(x_B, x_B)$. Since $\bar L > \bar V$, its slope is greater than 1.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 5. The feed condition $q$ and the q-line

    ### 5.1 Definition

    The feed changes the liquid flow at the feed stage. Define

    $$
    q \equiv \frac{\bar L - L}{F}
    $$

    that is, the moles of liquid that flow down the stripping section per mole of feed.
    Therefore

    $$
    \bar L = L + qF, \qquad \bar V = V - (1-q)\,F
    $$

    The second relation follows from the molar balance at the feed stage,
    $F + L + \bar V = \bar L + V$, after substituting $\bar L$.

    From an enthalpy balance, $q = \dfrac{H_V - H_F}{H_V - H_L}$
    ($H_V, H_L$ are saturated vapor and liquid enthalpies at the feed stage, $H_F$ the feed enthalpy).

    | Feed state | $q$ |
    |---|---|
    | Subcooled liquid | $q > 1$ |
    | Saturated liquid (bubble point) | $q = 1$ |
    | Partially vaporized | $0 < q < 1$ ($q$ = liquid fraction) |
    | Saturated vapor (dew point) | $q = 0$ |
    | Superheated vapor | $q < 0$ |

    ### 5.2 Derivation of the q-line

    At the feed stage both operating lines hold at the same point $(x, y)$. Write them in
    the "flow" form from Section 4:

    $$
    V y = L x + D x_D \qquad\text{(rectifying)}
    $$

    $$
    \bar V y = \bar L x - B x_B \qquad\text{(stripping)}
    $$

    Subtract the second from the first:

    $$
    (V - \bar V)\, y = (L - \bar L)\, x + D x_D + B x_B
    $$

    Now substitute $V - \bar V = (1-q)F$, $L - \bar L = -qF$ and $D x_D + B x_B = F z_F$.
    Divide by $F$:

    $$
    (1-q)\,y = -q\,x + z_F
    $$

    $$
    \boxed{\,q\,x + (1-q)\,y = z_F\,}
    \qquad\Longleftrightarrow\qquad
    y = \frac{q}{q-1}\,x - \frac{z_F}{q-1}
    $$

    **Checks:** at $x = z_F$ the left side is $q z_F + (1-q) y$, which equals $z_F$ when
    $y = z_F$, so the q-line passes through $(z_F, z_F)$. Its slope $q/(q-1)$ gives:

    - $q = 1$: **vertical** line $x = z_F$
    - $q = 0$: **horizontal** line $y = z_F$
    - $0 < q < 1$: negative slope
    - $q > 1$ or $q < 0$: positive slope

    We keep the form $q x + (1-q) y = z_F$ in the code because it has no division by zero at $q = 1$.

    ## 6. Locating the operating lines

    ### 6.1 Intersection point

    The two operating lines and the q-line meet at a single point $(x_i, y_i)$. Use the
    rectifying line $y = a x + b$ with $a = R/(R+1)$, $b = x_D/(R+1)$ together with the
    q-line, as a $2\times2$ linear system:

    $$
    \begin{aligned}
    a\,x - y &= -b \\
    q\,x + (1-q)\,y &= z_F
    \end{aligned}
    \qquad\Longleftrightarrow\qquad
    \begin{bmatrix} a & -1 \\ q & 1-q \end{bmatrix}
    \begin{bmatrix} x_i \\ y_i \end{bmatrix}
    =
    \begin{bmatrix} -b \\ z_F \end{bmatrix}
    $$

    ### 6.2 Stripping line

    Because the stripping line passes through $(x_B, x_B)$ and $(x_i, y_i)$, its slope is

    $$
    s = \frac{y_i - x_B}{x_i - x_B}
    \qquad\text{and}\qquad
    y = x_B + s\,(x - x_B)
    $$

    This slope equals $\bar L/\bar V$, which is a useful consistency check
    (the notebook prints both).

    ### 6.3 Flows

    $$
    L = R\,D, \quad V = (R+1)\,D, \quad \bar L = L + qF, \quad \bar V = \bar L - B
    $$
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 7. Stepping off stages, and limits

    ### 7.1 Stage stepping

    An **equilibrium stage** produces a vapor $y_n$ and a liquid $x_n$ in equilibrium.
    Two **passing streams** (vapor from stage $n+1$ and liquid from stage $n$) obey the operating line.
    On the diagram this is a staircase:

    1. Start at the top at $(x_D, x_D)$ and set $y_1 = x_D$ (total condenser).
    2. **Horizontal** move to the equilibrium curve: $x_n = \dfrac{y_n}{\alpha - (\alpha-1)y_n}$.
    3. **Vertical** move to the operating line to get $y_{n+1}$. Use the rectifying line
       while $x_n > x_i$. The first stage with $x_n \le x_i$ is the **feed stage**, and the
       stripping line is used from there on.
    4. Repeat until $x_n \le x_B$.

    Switching at the intersection gives the largest step possible, and therefore the
    **fewest** stages. This is the optimal feed location.

    The last stage is usually only partly needed. A fractional stage count is

    $$
    N = (n - 1) + \frac{x_{n-1} - x_B}{x_{n-1} - x_n}
    $$

    The count $N$ includes the reboiler (a partial reboiler is an equilibrium stage), so the column
    needs about $N - 1$ trays for a total condenser.

    ### 7.2 Minimum reflux

    As $R$ decreases, the intersection $(x_i, y_i)$ slides along the q-line toward the
    equilibrium curve. At the limit the operating lines touch the curve at the
    **pinch point** $(x_p, y_p)$, and an infinite number of stages is needed.

    - Pinch point: the q-line meets the equilibrium curve. For $q = 1$, $x_p = z_F$.
      Otherwise solve
      $\dfrac{\alpha x}{1+(\alpha-1)x} = \dfrac{q\,x - z_F}{q-1}$ (a quadratic; solved numerically here).
    - The rectifying line passes through $(x_D, x_D)$ and $(x_p, y_p)$, so its slope is

    $$
    \frac{R_{min}}{R_{min}+1} = \frac{x_D - y_p}{x_D - x_p}
    $$

    Let $s = \frac{x_D - y_p}{x_D - x_p}$. Then $R_{min} = s/(1-s)$, and since
    $1 - s = \frac{y_p - x_p}{x_D - x_p}$,

    $$
    \boxed{\,R_{min} = \frac{x_D - y_p}{y_p - x_p}\,}
    $$

    ### 7.3 Minimum stages (total reflux), Fenske equation

    At $R \to \infty$ both operating lines collapse onto $y = x$, so $y_{n+1} = x_n$. Each stage
    obeys the equilibrium relation in the form

    $$
    \frac{y_n}{1-y_n} = \alpha\,\frac{x_n}{1-x_n}
    $$

    Substituting $y_{n+1} = x_n$ repeatedly from the top ($y_1 = x_D$) to the bottom ($x_N = x_B$)
    multiplies by $\alpha$ once per stage:

    $$
    \frac{x_D}{1-x_D} = \alpha^{N_{min}}\,\frac{x_B}{1-x_B}
    \quad\Longrightarrow\quad
    \boxed{\,N_{min} = \frac{\ln\!\left[\dfrac{x_D}{1-x_D}\cdot\dfrac{1-x_B}{x_B}\right]}{\ln\alpha}\,}
    $$

    Together, $N_{min}$ and $R_{min}$ bracket a real design: a typical choice is $R = 1.05$ to $1.5\,R_{min}$.
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 8. Implementation

    Each function below is one equation from the derivation. Only the stage stepping uses a
    short loop, because each step depends on the previous one.
    """
    )
    return


@app.cell
def _(brentq, np):
    def y_eq(x, alpha):
        """Equilibrium curve y(x), constant relative volatility."""
        return alpha * x / (1 + (alpha - 1) * x)

    def x_eq(y, alpha):
        """Inverse equilibrium x(y): the horizontal step."""
        return y / (alpha - (alpha - 1) * y)

    def feed_intersection(zF, q, R, xD):
        """Rectifying line and q-line intersection as a 2x2 linear system (Section 6.1)."""
        a = R / (R + 1)
        b = xD / (R + 1)
        A = np.array([[a, -1.0],
                      [q, 1.0 - q]])
        rhs = np.array([-b, zF])
        xi, yi = np.linalg.solve(A, rhs)
        return xi, yi

    def min_reflux(zF, q, xD, alpha):
        """Pinch point where the q-line meets the equilibrium curve (Section 7.2)."""
        if q == 1:
            xp = zF
        else:
            xp = brentq(lambda x: y_eq(x, alpha) - (q * x - zF) / (q - 1), 1e-12, 1 - 1e-12)
        yp = y_eq(xp, alpha)
        return (xD - yp) / (yp - xp), xp, yp

    def min_stages(xD, xB, alpha):
        """Fenske equation (Section 7.3)."""
        return np.log(xD / (1 - xD) * (1 - xB) / xB) / np.log(alpha)

    def step_stages(zF, xD, xB, q, R, alpha, max_stages=200):
        """Step off stages from the top (Section 7.1).

        Returns the staircase vertices, fractional stage count, feed stage and intersection.
        """
        xi, yi = feed_intersection(zF, q, R, xD)
        slope_strip = (yi - xB) / (xi - xB)

        xs, ys = [xD], [xD]
        x_stage = []
        feed_stage = None
        y = xD
        N = np.inf
        for n in range(1, max_stages + 1):
            x = x_eq(y, alpha)                       # horizontal step
            x_stage.append(x)
            xs.append(x)
            ys.append(y)
            if feed_stage is None and x <= xi:
                feed_stage = n                       # switch operating line
            if x <= xB:
                x_prev = x_stage[-2] if n > 1 else xD
                N = (n - 1) + (x_prev - xB) / (x_prev - x)
                break
            if feed_stage is None:                   # vertical step
                y = R / (R + 1) * x + xD / (R + 1)
            else:
                y = xB + slope_strip * (x - xB)
            xs.append(x)
            ys.append(y)
        return np.array(xs), np.array(ys), N, feed_stage, xi, yi, slope_strip

    return (
        feed_intersection,
        min_reflux,
        min_stages,
        step_stages,
        x_eq,
        y_eq,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 9. Interactive explorer

    The reflux is set as a multiple of $R_{min}$, so the design is always feasible.
    """
    )
    return


@app.cell
def _(mo):
    alpha = mo.ui.slider(1.1, 6.0, step=0.1, value=2.5, label=r"$\alpha$", show_value=True)
    zF = mo.ui.slider(0.05, 0.95, step=0.05, value=0.5, label=r"$z_F$", show_value=True)
    xD = mo.ui.slider(0.60, 0.999, step=0.005, value=0.95, label=r"$x_D$", show_value=True)
    xB = mo.ui.slider(0.001, 0.40, step=0.005, value=0.05, label=r"$x_B$", show_value=True)
    q = mo.ui.slider(-0.5, 1.5, step=0.1, value=1.0, label=r"$q$", show_value=True)
    R_ratio = mo.ui.slider(1.02, 10.0, step=0.02, value=1.5, label=r"$R/R_{min}$", show_value=True)

    mo.hstack(
        [mo.vstack([alpha, zF, xD]), mo.vstack([xB, q, R_ratio])],
        justify="start",
        gap=2,
    )
    return R_ratio, alpha, q, xB, xD, zF


@app.cell
def _(R_ratio, alpha, min_reflux, min_stages, mo, q, step_stages, xB, xD, zF):
    mo.stop(
        not (xB.value < zF.value < xD.value),
        mo.md(r"Need $x_B < z_F < x_D$.").callout(kind="danger"),
    )

    _a, _zF, _xD, _xB, _q = alpha.value, zF.value, xD.value, xB.value, q.value

    Rmin, xp, yp = min_reflux(_zF, _q, _xD, _a)
    R = R_ratio.value * Rmin
    Nmin = min_stages(_xD, _xB, _a)
    xs, ys, N, feed_stage, xi, yi, slope_strip = step_stages(_zF, _xD, _xB, _q, R, _a)

    # Flows per 100 mol of feed
    F = 100.0
    D = F * (_zF - _xB) / (_xD - _xB)
    B = F - D
    L = R * D
    V = (R + 1) * D
    Lbar = L + _q * F
    Vbar = Lbar - B
    return B, D, F, L, Lbar, N, Nmin, R, Rmin, V, Vbar, feed_stage, slope_strip, xi, xp, xs, yi, yp, ys


@app.cell
def _(N, Nmin, R, Rmin, alpha, feed_stage, mo, np, plt, q, xB, xD, xi, xs, y_eq, yi, ys, zF):
    _x = np.linspace(0, 1, 400)
    _fig, _ax = plt.subplots(figsize=(6.5, 6.5))
    _ax.plot(_x, _x, "k--", lw=1)
    _ax.plot(_x, y_eq(_x, alpha.value), "b", lw=2, label="Equilibrium")
    _ax.plot([xD.value, xi], [xD.value, yi], "r", label="Rectifying line")
    _ax.plot([xB.value, xi], [xB.value, yi], "g", label="Stripping line")
    _ax.plot([zF.value, xi], [zF.value, yi], "m", label="q-line")
    _ax.plot(xs, ys, "k", lw=1.2)
    for _i, (_xv, _yv) in enumerate(zip(xs[1::2], ys[1::2]), start=1):
        _ax.annotate(str(_i), (_xv, _yv), textcoords="offset points", xytext=(-8, 4), fontsize=8)
    _ax.plot([xB.value, zF.value, xD.value], [xB.value, zF.value, xD.value], "ko", ms=4)
    _ax.set(
        xlim=(0, 1),
        ylim=(0, 1),
        xlabel="x (liquid mole fraction of A)",
        ylabel="y (vapor mole fraction of A)",
        title=f"N = {N:.1f} stages, feed on stage {feed_stage}, R = {R:.2f} "
              f"(Rmin = {Rmin:.2f}, Nmin = {Nmin:.1f})",
        aspect="equal",
    )
    _ax.legend(loc="lower right")
    _fig.tight_layout()
    _fig
    return


@app.cell
def _(B, D, F, L, Lbar, N, Nmin, R, Rmin, V, Vbar, feed_stage, mo, slope_strip, xi, xp, yi, yp):
    mo.md(
        f"""
    | Result | Value |
    |---|---|
    | Stages $N$ (incl. reboiler) | {N:.2f} |
    | Optimal feed stage (from top) | {feed_stage} |
    | Minimum stages (Fenske) | {Nmin:.2f} |
    | Minimum reflux | {Rmin:.3f} |
    | Operating reflux | {R:.3f} |
    | Pinch point $(x_p, y_p)$ | ({xp:.3f}, {yp:.3f}) |
    | Intersection $(x_i, y_i)$ | ({xi:.3f}, {yi:.3f}) |
    | Flows per {F:.0f} mol feed: $D$, $B$ | {D:.1f}, {B:.1f} |
    | $L$, $V$ | {L:.1f}, {V:.1f} |
    | $\\bar L$, $\\bar V$ | {Lbar:.1f}, {Vbar:.1f} |
    | Stripping slope: geometric vs $\\bar L/\\bar V$ | {slope_strip:.4f} vs {Lbar / Vbar:.4f} |
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## 10. Exercises

    Use the sliders to answer these, then check the reasoning against the derivation.

    1. Reduce $R/R_{min}$ toward 1. What happens to $N$, and why does the pinch region matter?
    2. Set $R/R_{min}$ very large. Does $N$ approach $N_{min}$? Compare with the Fenske result.
    3. Move $q$ from 1.5 to $-0.5$. How does the intersection move along the q-line, and what happens to the feed stage?
    4. Increase $\alpha$ from 1.5 to 5. Which of $N$, $R_{min}$ and $N_{min}$ change most?
    5. Set $q = 0.5$ and verify by hand that the q-line goes through $(z_F, z_F)$ with slope $q/(q-1) = -1$.
    6. Verify the table row "Stripping slope": explain why $\bar L/\bar V$ equals the geometric slope.
    7. **Code exercise:** replace `y_eq` and `x_eq` with tabulated equilibrium data using
       `scipy.interpolate.PchipInterpolator` (and its inverse) to handle a non-ideal system.
    8. **Code exercise:** add a Murphree vapor efficiency by stepping only a fraction
       $E_{MV}$ of the way to the equilibrium curve on each stage.

    ## Limitations

    Constant relative volatility, constant molal overflow, total condenser, a single feed,
    no side streams and no heat losses. Azeotropic or strongly non-ideal systems need
    real equilibrium data (see exercise 7).

    To share it as a static in-browser (WebAssembly) page:
    `marimo export html-wasm mccabe_thiele.py -o dist --mode run`
    """
    )
    return


if __name__ == "__main__":
    app.run()
