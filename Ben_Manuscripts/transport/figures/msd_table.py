#!/usr/bin/env python

import sqlite3 as sql

connection = sql.connect("../../../timeseries/msd.db")
crsr = connection.cursor()
query = "SELECT name, sigma, alpha, hurst, python_MSD, python_MSD_CI_lower, python_MSD_CI_upper from msd ORDER BY python_MSD DESC"
output = crsr.execute(query).fetchall()

# dictionary of residue names and their corresponding chemical names as they will be dispalyed in the table
res_to_name = {}
res_to_name["ACN"] = "Acetamide"
res_to_name["ACH"] = "Acetic Acid"
res_to_name["ATO"] = "Acetone"
res_to_name["BUT"] = "Butanol"
res_to_name["DMF"] = "Dimethyl Formamide"
res_to_name["DMS"] = "Dimethyl Sulfoxide"
res_to_name["DMP"] = "2,3-Dimercapto-1-propanol"
res_to_name["EAC"] = "Ethyl Acetate"
res_to_name["ETH"] = "Ethanol"
res_to_name["GCL"] = "Ethylene Glycol"
res_to_name["GLY"] = "Glycerol"
res_to_name["MET"] = "Methanol"
res_to_name["PR"] = "Propanol"
res_to_name["PCB"] = "Propylene Carbonate"
res_to_name["PG"] = "Propylene Glycol"
res_to_name["RIB"] = "Ribose"
res_to_name["SOH"] = "Mercaptoethanol"
res_to_name["TET"] = "Tetrose"
res_to_name["THF"] = "Tetrahydrofuran"
res_to_name["URE"] = "Urea"

table = ""
table += r"\begin{table}[h]" + "\n" + r"  \centering" + "\n" + r"  \begin{tabular}{ccccc}" + "\n" + r"  \toprule" + "\n"
table += r"  System & $\sigma$ ($nm$) & $\alpha$ & $H$ & Simulated MSD ($nm^2$)\\" + "\n" + r"  \midrule" + "\n"

for i in output:
	table += r"  %s & %.2f & %.2f & %.2f & %.2f [%.2f, %.2f] \\" % (res_to_name[i[0]], i[1], i[2], i[3], i[4], i[5], i[6]) + "\n"


table += r"  \bottomrule" + "\n" + r"  \end{tabular}"
print(table)
