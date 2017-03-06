#!/usr/bin/env python
import pexpect

prompts = [">>> ",""]
process = pexpect.spawn("python")
process.expect(prompts)
process.sendline("print 'hello'" )
index = process.expect(prompts)
print 1,process.before, process.after
process.sendline("")
index = process.expect(prompts)
print 2,process.before, process.after


