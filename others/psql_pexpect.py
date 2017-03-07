#!/usr/bin/env python
import pexpect

process = pexpect.spawn("psql -U contrax ddtcontrax -p 5433")
process.expect(".*=# ")
process.sendline("\d")
result = ""
i=process.expect([":", "(END).*", ".*=# .*"])		#Prompt que indica que esta esperando por nueva orden; pueden ser varios, y la i indica cual ha  leido
while (i == 0):
	print process.before										# En process.before se almacena el resultado antes del prompt
	process.send(" ")
	i=process.expect([":", "(END).*", ".*=# .*"])

process.send("q")
process.expect(".*=# ")
process.sendline("select id,code from \"KER_thirdpart\" limit 100;")
i=process.expect([":", "(END).*", ".*=# .*"])
while (i == 0):
	print process.before
	process.send(" ")
	i=process.expect([":", "(END).*", ".*=# .*"])

process.send("q")
