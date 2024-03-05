# Read the file containing instance names
$instanceNames = Get-Content -Path "instance_names.txt"

# Create a new CSV file for output
$outputFilePath = "output.csv"
"Instance,Cost,Time" | Out-File -FilePath $outputFilePath -Encoding ASCII

# Iterate through each instance name
foreach ($instanceName in $instanceNames) {
    # Construct the argument with instance name
    $argument = ".\instances\$instanceName.tsp"

    Write-Host "Running instance: $instanceName"
    Write-Host "File path: $argument"

    # Run main.exe with the argument
    $output = ./main.exe $argument

    # Split the output into lines
    $lines = $output -split "`n"

    # Initialize cost and cpuExecutionTime as "N/A"
    $cost = "N/A"
    $cpuExecutionTime = "N/A"

    # Loop through each line
    foreach ($line in $lines) {
        # Extract the cost
        if ($line -match 'Average s cost: ([\d.]+)') {
            $cost = $matches[1]
        }

        # Extract the average CPU execution time
        if ($line -match 'Average CPU execution time: ([\d.]+) s') {
            $cpuExecutionTime = $matches[1]
        }
    }

    # Output instance name, cost, and CPU execution time to the CSV file
    "$instanceName,$cost,$cpuExecutionTime" | Out-File -FilePath $outputFilePath -Encoding ASCII -Append
}