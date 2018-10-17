import subprocess
bashCommand = "sage test_for_sage.sage"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
print output
