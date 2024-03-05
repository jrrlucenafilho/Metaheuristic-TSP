# Read the file containing instance names
$instanceNames = Get-Content -Path "instance_names.txt"

# Create a new CSV file for output
$outputFilePath = "output.csv"
"Instance,Cost,AvgTime,BuildSolution,BestImprovementSwap,OrOpt,OrOpt2,OrOpt3,2-Opt,Disturbance" | Out-File -FilePath $outputFilePath -Encoding ASCII

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

    # Initialize variables
    $cost = "N/A"
    $cpuExecutionTime = "N/A"
    $buildSolution = "N/A"
    $bestImprovementSwap = "N/A"
    $orOpt = "N/A"
    $orOpt2 = "N/A"
    $orOpt3 = "N/A"
    $twoOpt = "N/A"
    $disturbance = "N/A"

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

        # Extract the times for each neighborhood structure
        if ($line -match 'BuildSolution time: ([\d.]+) s') {
            $buildSolution = $matches[1]
        }

        if ($line -match 'BestImprovementSwap time: ([\d.]+) s') {
            $bestImprovementSwap = $matches[1]
        }

        if ($line -match 'OrOpt time: ([\d.]+) s') {
            $orOpt = $matches[1]
        }

        if ($line -match 'OrOpt2 time: ([\d.]+) s') {
            $orOpt2 = $matches[1]
        }

        if ($line -match 'OrOpt3 time: ([\d.]+) s') {
            $orOpt3 = $matches[1]
        }

        if ($line -match '2-Opt time: ([\d.]+) s') {
            $twoOpt = $matches[1]
        }

        if ($line -match 'Disturbance time: ([\d.]+) s') {
            $disturbance = $matches[1]
        }
    }

    # Output instance name, cost, CPU execution time, and times for each neighborhood structure to the CSV file
    "$instanceName,$cost,$cpuExecutionTime,$buildSolution,$bestImprovementSwap,$orOpt,$orOpt2,$orOpt3,$twoOpt,$disturbance" | Out-File -FilePath $outputFilePath -Encoding ASCII -Append
}