import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas, re, glob
import pynamd

def Hill(pH, n, pKa):
    S = 1 / ( (10 ** (n*(pKa - pH))) + 1)
    return S

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

root = "examples"

configfiles = [f"{root}/Radak.json"]
cphlogs = natural_sort([x.replace("\\", "/") for x in glob.glob(f"{root}/pH*/*cphlog")])

S_over_time = {}
for cphlog in cphlogs:

    if "_" in cphlog.split("pH")[1].split("/")[0]:
        run = int(cphlog.split("pH")[1].split("/")[0].split("_")[1])
    else:
        run = 0
        
    pynamd_parser = pynamd.cphlog.TitratableSystemSet.from_cphlogs([cphlog], configfiles)
    pH = float(pynamd_parser.pHs[0])

    states = pynamd_parser.macro_occupancies() #1 = proton occupies space
    # so we reverse it to the S value which is 1 = deprotonated
    states -= 1
    states = np.abs(states)
    groups = pynamd_parser.segresids()
    if run == 0:
        #S_over_time[pH] = pandas.DataFrame(states, np.arange(1, states.shape[0]+1), columns=groups)
        S_over_time[pH] = states
    else:
        S_over_time[pH] =  np.vstack((S_over_time[pH], states))


X = []
Y = []
for pH in S_over_time.keys():
    #Titration = pandas.DataFrame(S_over_time[pH], columns=groups)
    S = S_over_time[pH][-1000].mean() # take last 1000 logs for the mean
    X.append(pH)
    Y.append(S.mean())
    print(f"pH {pH}:", S)

plt.figure(figsize=(5,3))

plt.plot(np.linspace(min(X), max(X), 100), Hill(np.linspace(min(X), max(X), 100), 1, 6.0), label="Ref His")

plt.scatter(X, Y, color="orange", marker="1")
popt, _pcov = curve_fit(Hill, X, Y, bounds=[(0,3.0), (1, 14)])
n = popt[0].round(2)
pKa = popt[1].round(2)

plt.plot(np.linspace(min(X), max(X), 100), Hill(np.linspace(min(X), max(X), 100), n, pKa), color="orange",label="Measured His")

plt.text(6.5, 0.5, f"Apparent pKa: {pKa}")
plt.text(6.51, 0.4, f"Cooperativity: {n}")

plt.ylabel("S")
plt.xlabel("pH")
plt.legend()
print(f"pKa = {pKa} n = {n}")

plt.title("Histidine apparent pKa in a\nHYF tripeptide crystal at 50% water")
plt.tight_layout()
fig1 = plt.gcf()
plt.show()
plt.draw()
fig1.savefig(f"CpHMD.png", dpi=400, bbox_inches='tight')