#!/usr/bin/env python

arq = open("Bonamente_table.dat", "r")

outPut = open("fractionBarionsSZ.dat", "w")

def decompoeDado(dado):
    c = dado.replace("$", "")
    c = c.replace(" ", "")
    c = c.replace("\\\\", "")
    c = c.replace("{", "")
    c =c.replace("}", "")
    c = c.replace("\\pm", ",")
    c = c.replace("^", "")
    c = c.replace("_", ",")

    return c


for li in arq.readlines():
    b = li.split("&")
    c = b[8]
    d = b[10]
    e = b[1].replace(" ", "")+","+decompoeDado(c)+","+decompoeDado(d)
    print(e)
    outPut.write(e)
 
outPut.close()
