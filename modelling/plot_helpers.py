from cycler import cycler


def custom_aptavita_colors():
    return ["#9B0138", "#057D54", "#FFCE3A", "#EDB4D7", "#4FD590", "#4D94EF", "#EFA54D", "#313175",
            "#667817", "#6FC0A8", "#D6682A", "#F3758A", "#755F26", "#A74D36", "#C88E99", "#A79536"]


def custom_aptavita_color_cycler():
    colors = custom_aptavita_colors()
    return cycler(color=colors)


def senstivity_analysis_factor_names():
    return ["$k_{\mathrm{ts}}$", "$k_{\mathrm{tl}}$", "$k_{\mathrm{mat}}$", "$k_{\mathrm{cat}}$", "$K_{\mathrm{s}}$", "$kc_{\mathrm{s}}$", "$K_{\mathrm{l}}$", "$K_{\mathrm{TlR}}$", "$K_{\mathrm{m}}$", "$\delta_{\mathrm{mRNA}}$", "$\delta_{\mathrm{TlR}}$", "$k_{\mathrm{on}}$", "$k_{\mathrm{off}}$", "$k_{\mathrm{c}}$", "$\mathrm{DNA}$", "$\mathrm{Vit}_\mathrm{tot}$"]


def add_value_labels(ax, spacing=5):
    """Add labels to the end of each bar in a bar chart.

    Arguments
    ---------
        ax: matplotlib.axes.Axes
            The matplotlib object containing the axes of the plot to annotate.
        spacing: int
            The distance between the labels and the bars. (default 5)

    Notes
    -----
    Function copied and adapted from: https://stackoverflow.com/questions/28931224/adding-value-labels-on-a-matplotlib-bar-chart
    """

    # For each bar: Place a label
    for rect in ax.patches:
        # Get X and Y placement of label from rect.
        y_value = rect.get_height()
        x_value = rect.get_x() + rect.get_width() / 2

        # Number of points between bar and label. Change to your liking.
        space = spacing
        # Vertical alignment for positive values
        va = 'bottom'

        # If value of bar is negative: Place label below bar
        if y_value < 0:
            # Invert space to place label below
            space *= -1
            # Vertically align label at top
            va = 'top'

        # Use Y value as label and format number with one decimal place
        label = "{:.1f}".format(y_value)

        # Create annotation
        ax.annotate(
            label,                      # Use `label` as label
            (x_value, y_value),         # Place label at end of the bar
            xytext=(0, space),          # Vertically shift label by `space`
            textcoords="offset points",  # Interpret `xytext` as offset in points
            ha='center',                # Horizontally center label
            va=va)                      # Vertically align label differently for
        # positive and negative values.
