import os
import argparse
from rich.table import Table
from rich.console import Console

# parse arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("--analysis", type=str, help="name of the analysis")
parser.add_argument("--config", type=str, help="name of the config")
parser.add_argument("--samplelist", type=str, help="name of the samplelist")
parser.add_argument("--tag", type=str, help="name of the tag")

args = parser.parse_args()

# build the law command to be run to get the status of the production
law_cmd = f"law run ProduceSamples --analysis {args.analysis} --config {args.config} --sample-list {args.samplelist} --production-tag {args.tag} --workers 1 --print-status 1"
# now run the law command using subprocess and store the output in a variable
output = os.popen(law_cmd).read()
data = {}
# now parse output line by line
lines = output.split("\n")
for i, line in enumerate(lines):
    parsing_sample = False
    # if CROWNRun is in the line, we get the status of a new sample
    if not parsing_sample:
        if "> CROWNRun" in line:
            parsing_sample = True
            # find out the samplename
            samplename = line.split("nick=")[1].split(",")[0]
            data[samplename] = {}
        else:
            if "NestedSiblingFileCollection" in line:
                statusline = lines[i + 1]
                result = statusline[statusline.find("(") + 1 : statusline.find(")")]
                data[samplename]["done"] = int(result.split("/")[0])
                data[samplename]["total"] = int(result.split("/")[1])
                parsing_sample = False

# now print the results as a rich table with the percentual completion

table = Table(title="Sample Status")

table.add_column("Sample", justify="right")
table.add_column("Done", justify="right")
table.add_column("Total", justify="right")
table.add_column("Percent", justify="right")


for sample in data:
    done = data[sample]["done"]
    total = data[sample]["total"]
    percent = int(float(done) / float(total) * 100)
    table.add_row(sample, str(done), str(total), str(percent) + "%")
# add a total row at the end with the sum of all percentual completion
total_done = sum([data[sample]["done"] for sample in data])
total_total = sum([data[sample]["total"] for sample in data])
table.add_row(
    "Total",
    str(total_done),
    str(total_total),
    str(int(float(total_done) / float(total_total) * 100)) + "%",
)

console = Console()
console.print(table)
