import gffutils
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


def plot_coverage(input, gtf, output):
    db = gffutils.create_db(gtf, ':memory:')
    data = pd.read_table(input, header=None)

    line_height = 14
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(data[1], data[2], c="b", label="sense")
    ax.plot(data[1], data[5], c="r", label="antisense")
    ax.legend()
    lines = []
    fig.canvas.draw()
    fig.subplots_adjust(bottom=0.3)
    data_fig = ax.transData + fig.transFigure.inverted()
    line_height = fig.transFigure.inverted().transform([0, line_height])[1]
    for f in db.features_of_type("exon"):
        start = data_fig.transform([f.start, 0])[0]
        end = data_fig.transform([f.end, 0])[0]

        r = mpl.patches.Rectangle(
            (start, 0.3 - 2.5 * line_height),
            end - start,
            line_height * 0.4,
            fill=True,
            facecolor="#90a0b0",
            linewidth=1,
            zorder=5,
            edgecolor="white"
        )
        fig.add_artist(r)
        rx, ry = r.get_xy()
        cx = rx + r.get_width() / 2
        cy = ry + r.get_height() / 2

        ann = mpl.text.Text(
            cx, cy,
            f.attributes["gene_name"][0],
            color="black",
            ha="center",
            va="center",
            fontsize=10,
            zorder=10
        )
        fig.add_artist(ann)

        rect_bb = r.get_window_extent()
        rect_bb = rect_bb.transformed(fig.transFigure.inverted())
        text_bb = ann.get_window_extent()
        text_bb = text_bb.transformed(fig.transFigure.inverted())
        if len(lines) == 0:
            lines.append((rect_bb, text_bb))
        else:
            found = False
            for i, prev_bb in enumerate(lines):
                if not text_bb.overlaps(prev_bb[0]) and not text_bb.overlaps(prev_bb[1]):
                    lines[i] = rect_bb, text_bb
                    found = True
                    break
            if not found:
                i = len(lines)
                lines.append((rect_bb, text_bb))
            if i > 0:
                r.set_y(r.get_y() - line_height * i * 2)
            ann.set_position((cx, cy - line_height * i * 2))
        ann_pos = ann.get_position()
        ann.set_position((ann_pos[0], ann_pos[1] - 1 * line_height))
    fig.savefig(output)
