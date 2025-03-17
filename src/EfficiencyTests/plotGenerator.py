import matplotlib.pyplot as plt
import sys
import ast  # For safely evaluating strings as Python literals

# This Python script will read the .txt file given and make a plot
# The .txt files from the folder "results" must have a specific format:
# 1st line: Title of the plot
# 2nd line: Label of the first algorithm used
# 3rd line: Results written as an array
# 4th line: Label of the second algorithm used
# 5th line: Results written as an array
# etc...

def plotGenerator(file_path, timeOrMem, fixedVariable, linOrExp):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Title
    title = lines[0].strip()

    algorithms = []
    results = []

    i = 1
    while i < len(lines):
        algorithm_name = lines[i].strip()
        algorithms.append(algorithm_name)

        results_array = ast.literal_eval(lines[i + 1].strip())
        results.append(results_array)

        # Next algorithm
        i += 2

    plt.figure()

    # Determine x_values based on linOrExp
    if linOrExp == 'exp':
        x_values = [2**i for i in range(1, 16)]
    else:
        x_values = [1 + 200 * i for i in range(101)]

    for algo, res in zip(algorithms, results):
        plt.plot(x_values, res, label=algo, marker='.')

    # Second argument when calling the function should be "time" or "memory"
    if timeOrMem == "time":
        plt.ylabel('Time (ticks)')
    else:
        plt.ylabel('Memory (bytes)')

    if fixedVariable == '0':
        plt.xlabel('Degree')  # if fixedVariable == 0: the degree is changing and the coeff size is fixed
    else:
        plt.xlabel('Coefficient size (bits)')  # if fixedVariable == 1: the coeff size is changing and the degree is fixed

    plt.xscale("log")  # Logarithmic scale for better visualization
    plt.yscale("log")
    plt.title(title)
    plt.legend()
    plt.grid(True)

    save_file_path = file_path.rsplit('.', 1)[0] + '.png'
    plt.savefig(save_file_path)

# Main script logic
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python plot_results.py <file_path> <'time' or 'memory'> <0 or 1> <'lin' or 'exp'>")
        sys.exit(1)

    file_path = sys.argv[1]
    timeOrMem = sys.argv[2]
    fixedVariable = sys.argv[3]
    linOrExp = sys.argv[4]

    plotGenerator(file_path, timeOrMem, fixedVariable, linOrExp)
