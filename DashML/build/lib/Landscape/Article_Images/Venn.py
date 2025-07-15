import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

venn3(subsets=(100, 100,40, 100,37, 40,20), set_labels = ('Nanopore', 'ShapeMap', 'ShapeCE'))
# A, B, AB, C, AC, BC,ABC
plt.savefig(fname="venn_diagram.png", dpi=600)
#plt.show)
