from NETS import NETS
from NETS import Tuple

amount = {}
amount[1] = 3
amount[7] = 5
amount[8] = 5
id = 0
windowSize = 3
slideSize = 1
Window = {}
nets = NETS(dim=1, subDim=1, R=2.0, k=6, S=slideSize, W=windowSize, nW=1, maxValues=[8], minValues=[1])
for j in range(windowSize):
    for val in amount:
        for i in range(amount[val]):
            if j not in Window:
                Window[j] = []
            data = Tuple(id, j, [val])
            id += 1
            Window[j].append(data)
    nets.calcNetChange(Window[j], j)
    nets.findOutlier(j)
    outliers = nets.returnOutliers()
    print(outliers)