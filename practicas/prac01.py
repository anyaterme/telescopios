from math import *
import numpy as np

reflection = {}

reflection["lente"] = [1289,1300,1318,1325]
#reflection["mirrow_s"] = [38,38,38,38]
reflection["glass"] = [1160,1090,1154,1157]
reflection["oro"] = [556,544,597,590]
reflection["al_xl"] = [890,884,884,891,895]
reflection["al_sm"] = [690,673,672,672,671]
reflection["opaco"] = [73,73,73,73,73]
reflection["air"] = [1390,1383,1380,1390,1385,1390,1390,1390,1390,1390,1386,1390,1400,1394,1398,1407,1390,1390,1390]


nolaser = 34

mean = {}
for key in reflection:
	mean[key] = (np.mean(reflection[key]) - nolaser) * 0.001

complementary = ["lente", "glass", "air"]

for key in mean:
	if key in complementary:
		print "%s refleja el %lf %%" % (key,((1 - (mean[key]/mean["air"]))*100.))
	else:
		print "%s refleja el %lf %%" % (key,(((mean[key]/mean["air"]))*100.))
