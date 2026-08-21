plugins {
    base
}

group = "pbgz.ci"
version = "2.2.0"

val workspaceDir = file(layout.projectDirectory.dir("workspace"))
val logsDir = file(layout.projectDirectory.dir("logs"))
val timestamp = java.time.LocalDateTime.now().format(java.time.format.DateTimeFormatter.ofPattern("yyyyMMdd_HHmmss"))

repositories {
    mavenCentral()
}

tasks.register("cleanWorkspace", Delete::class) {
    description = "Clean workspace"
    group = "ci"
    delete(workspaceDir)
}

tasks.register("setupDirectories") {
    description = "Create necessary directories"
    group = "ci"
    
    doLast {
        workspaceDir.mkdirs()
        logsDir.mkdirs()
        println("[$timestamp] directory creation complete")
        println("workspace directory: ${workspaceDir.absolutePath}")
        println("log directory: ${logsDir.absolutePath}")
    }
}

tasks.register("downloadSource") {
    description = "Download source code from GitHub (with retry and integrity check)"
    group = "ci"
    
    dependsOn("setupDirectories")
    
    val githubRepo = project.findProperty("github.repo") as String? ?: ""
    val githubBranch = project.findProperty("github.branch") as String? ?: "v2.2.0"
    val maxAttempts = 3
    
    if (githubRepo.isEmpty()) {
        throw GradleException("Error: github.repo not configured, please set it in gradle.properties")
    }
    
    val sourceCodeDir = file("${workspaceDir}/pbgz")
    
    doLast {
        if (sourceCodeDir.exists()) {
            println("[$timestamp] source code directory already exists, cleaning up...")
            delete(sourceCodeDir)
        }
        
        var success = false
        var attempt = 1
        
        while (attempt <= maxAttempts && !success) {
            println("[$timestamp] === Attempting to download source code (attempt ${attempt}/${maxAttempts}) ===")
            
            try {
                // execute git clone
                val process = ProcessBuilder("git", "clone", "-b", githubBranch, githubRepo, sourceCodeDir.absolutePath)
                    .directory(project.projectDir)
                    .inheritIO()
                    .start()
                
                val exitCode = process.waitFor()
                
                if (exitCode == 0) {
                    println("[$timestamp] download complete")
                    
                    // check source code integrity
                    println("[$timestamp] checking source code integrity...")
                    
                    if (!sourceCodeDir.exists()) {
                        println("[$timestamp] ❌ Error: source code directory does not exist")
                        delete(sourceCodeDir)
                    } else if (!file("${sourceCodeDir}/.git").exists()) {
                        println("[$timestamp] ❌ Error: .git directory does not exist, download may be incomplete")
                        delete(sourceCodeDir)
                    } else if (!file("${sourceCodeDir}/README.md").exists() || 
                              !file("${sourceCodeDir}/CMakeLists.txt").exists() ||
                              !file("${sourceCodeDir}/build-release.sh").exists()) {
                        println("[$timestamp] ❌ Error: key files do not exist, download may be incomplete")
                        delete(sourceCodeDir)
                    } else {
                        // check git repository status
                        val gitCheck = ProcessBuilder("git", "rev-parse", "--git-dir")
                            .directory(sourceCodeDir)
                            .start()
                        
                        val gitExit = gitCheck.waitFor()
                        
                        if (gitExit == 0) {
                            // check file count (at least 100 files)
                            val fileCount = ProcessBuilder("find", ".", "-type", "f")
                                .directory(sourceCodeDir)
                                .start()
                            
                            val fileCountText = fileCount.inputStream.bufferedReader().use { it.readText() }
                            val fileCountLines = fileCountText.lines().size
                            
                            if (fileCountLines >= 100) {
                                success = true
                                println("[$timestamp] ✅ source code download and verify successful")
                                
                                // display version information
                                val gitLog = ProcessBuilder("git", "log", "-1", "--oneline", "--format=%H - %s (%an, %ar)")
                                    .directory(sourceCodeDir)
                                    .start()
                                
                                val versionInfo = gitLog.inputStream.bufferedReader().use { it.readText() }
                                println("[$timestamp] current version information: $versionInfo")
                            } else {
                                println("[$timestamp] ❌ Error: abnormal file count (${fileCountLines} files), download may be incomplete")
                                delete(sourceCodeDir)
                            }
                        } else {
                            println("[$timestamp] ❌ Error: git repository status abnormal")
                            delete(sourceCodeDir)
                        }
                    }
                } else {
                    println("[$timestamp] ❌ git clone failed (exit code: $exitCode)")
                    delete(sourceCodeDir)
                }
            } catch (e: Exception) {
                println("[$timestamp] ❌ Exception occurred during download: ${e.message}")
                delete(sourceCodeDir)
            }
            
            if (!success && attempt < maxAttempts) {
                val waitTime = attempt * 5
                println("[$timestamp] waiting ${waitTime} seconds before retry...")
                Thread.sleep(waitTime * 1000L)
            }
            
            attempt++
        }
        
        if (!success) {
            throw GradleException("Error: source code download still failed after ${maxAttempts} attempts. Please check network connection and configuration before retry.")
        }
        
        println("[$timestamp] source code download complete (branch: $githubBranch)")
        
        val logFile = file("${logsDir}/download_source_${timestamp}.log")
        logFile.writeText("source code download complete - branch: $githubBranch\ntime: $timestamp")
    }
}

tasks.register<Exec>("buildGcc") {
    description = "Run build.all.gcc.sh build script"
    group = "ci"
    
    dependsOn("downloadSource")
    
    val buildScriptGcc = project.findProperty("build.script.gcc") as String? ?: "3rd_party/build.all.gcc.sh"
    val sourceCodeDir = file("${workspaceDir}/pbgz")
    val buildScriptPath = file("${sourceCodeDir}/${buildScriptGcc}")
    
    doFirst {
        if (!buildScriptPath.exists()) {
            throw GradleException("Error: build script does not exist: ${buildScriptPath.absolutePath}")
        }
        if (!buildScriptPath.canExecute()) {
            println("[$timestamp] Warning: script not executable, adding execute permission...")
            buildScriptPath.setExecutable(true)
        }
        // Ensure log directory exists before opening log file stream
        logsDir.mkdirs()
    }
    
    commandLine = mutableListOf("bash", buildScriptPath.absolutePath)
    workingDir = sourceCodeDir
    
    val logFile = file("${logsDir}/build_gcc_${timestamp}.log")
    doFirst {
        standardOutput = java.io.FileOutputStream(logFile)
        errorOutput = System.err
    }
    
    doLast {
        println("[$timestamp] build.all.gcc.sh execution complete")
    }
}

tasks.register<Exec>("buildRelease") {
    description = "Run build-release.sh build script"
    group = "ci"
    
    dependsOn("buildGcc")
    
    val buildScriptRelease = project.findProperty("build.script.release") as String? ?: "build-release.sh"
    val sourceCodeDir = file("${workspaceDir}/pbgz")
    val buildScriptPath = file("${sourceCodeDir}/${buildScriptRelease}")
    
    doFirst {
        if (!buildScriptPath.exists()) {
            throw GradleException("Error: build script does not exist: ${buildScriptPath.absolutePath}")
        }
        if (!buildScriptPath.canExecute()) {
            println("[$timestamp] Warning: script not executable, adding execute permission...")
            buildScriptPath.setExecutable(true)
        }
        // Ensure log directory exists before opening log file stream
        logsDir.mkdirs()
    }
    
    commandLine = mutableListOf("bash", buildScriptPath.absolutePath)
    workingDir = sourceCodeDir
    
    val logFile = file("${logsDir}/build_release_${timestamp}.log")
    doFirst {
        standardOutput = java.io.FileOutputStream(logFile)
        errorOutput = System.err
    }
    
    doLast {
        println("[$timestamp] build-release.sh execution complete")
    }
}

tasks.register<Exec>("runTests") {
    description = "Run tests"
    group = "ci"
    
    dependsOn("buildRelease")
    
    val runTests = project.findProperty("test.enabled") as String? ?: "false"
    
    doFirst {
        if (runTests != "true") {
            println("[$timestamp] skipping tests (test.enabled=false)")
            throw StopExecutionException("tests disabled")
        }
    }
    
    val testScript = project.findProperty("test.script") as String? ?: ""
    val testCommand = project.findProperty("test.command") as String? ?: ""
    val sourceCodeDir = file("${workspaceDir}/pbgz")
    
    if (testScript.isNotEmpty()) {
        val testScriptPath = file("${sourceCodeDir}/${testScript}")
        doFirst {
            if (!testScriptPath.exists()) {
                throw GradleException("Error: test script does not exist: ${testScriptPath.absolutePath}")
            }
            
            // Check if it's a Python script
            if (testScript.endsWith(".py")) {
                // Check Python availability
                val pythonCheck = ProcessBuilder("python3", "--version")
                    .inheritIO()
                    .start()
                if (pythonCheck.waitFor() != 0) {
                    throw GradleException("Error: python3 not available for running Python test scripts")
                }
                println("[$timestamp] executing Python test script: ${testScriptPath.name}")
            } else {
                if (!testScriptPath.canExecute()) {
                    println("[$timestamp] Warning: test script not executable, adding execute permission...")
                    testScriptPath.setExecutable(true)
                }
            }
        }
        
        if (testScript.endsWith(".py")) {
            commandLine = mutableListOf("python3", testScriptPath.absolutePath)
        } else {
            commandLine = mutableListOf("bash", testScriptPath.absolutePath)
        }
        workingDir = sourceCodeDir
    } else if (testCommand.isNotEmpty()) {
        val commandParts = testCommand.split(" ".toRegex())
        commandLine = commandParts.toMutableList()
        workingDir = sourceCodeDir
    } else {
        doFirst {
            println("[$timestamp] Warning: test script or command not configured")
            throw StopExecutionException("tests not configured")
        }
    }
    
    val logFile = file("${logsDir}/test_${timestamp}.log")
    doFirst {
        logsDir.mkdirs()
        standardOutput = java.io.FileOutputStream(logFile)
        errorOutput = System.err
    }
    
    doLast {
        println("[$timestamp] test execution complete")
        println("[$timestamp] test log saved to: ${logFile.absolutePath}")
    }
}

tasks.register("ci") {
    description = "Run complete CI process: download, build, test"
    group = "ci"
    
    val notifyEnabled = project.findProperty("notify.enabled") as String? ?: "true"
    val notifyScript = project.findProperty("notify.script") as String? ?: "./notify.sh"
    val ciLogFile = file("${logsDir}/ci_${timestamp}.log")
    
    doFirst {
        if (notifyEnabled == "true" && File(notifyScript).exists()) {
            println("[$timestamp] sending build start notification...")
            exec {
                commandLine = mutableListOf("bash", notifyScript, "start")
                workingDir = project.projectDir
            }
        }
    }
    
    dependsOn("setupDirectories")
    dependsOn("downloadSource")
    dependsOn("buildGcc")
    dependsOn("buildRelease")
    
    val runTests = project.findProperty("test.enabled") as String? ?: "false"
    if (runTests == "true") {
        dependsOn("runTests")
    }
    
    doLast {
        println("=== CI execution summary ===")
        println("working directory: ${workspaceDir.absolutePath}")
        println("log directory: ${logsDir.absolutePath}")
        println("execution time: ${java.time.LocalDateTime.now()}")
        println("status: success")
        
        val latestLog = if (ciLogFile.exists()) ciLogFile.absolutePath else "log file not generated"
        
        if (notifyEnabled == "true" && File(notifyScript).exists()) {
            println("[$timestamp] sending build success notification...")
            exec {
                commandLine = mutableListOf("bash", notifyScript, "success", latestLog)
                workingDir = project.projectDir
            }
        }
    }
}

tasks.register("testNotification") {
    description = "Test notification system"
    group = "help"
    
    val notifyScript = project.findProperty("notify.script") as String? ?: "./notify.sh"
    
    doLast {
        if (!File(notifyScript).exists()) {
            throw GradleException("notification script does not exist: $notifyScript")
        }
        
        println("[$timestamp] testing notification system...")
        exec {
            commandLine = mutableListOf("bash", notifyScript, "test")
            workingDir = project.projectDir
        }
    }
}

tasks.register("showConfig") {
    description = "Display current configuration"
    group = "help"
    
    doLast {
        println("=== current CI configuration ===")
        println("GitHub repository: ${project.findProperty("github.repo") ?: "not configured"}")
        println("branch: ${project.findProperty("github.branch") ?: "not configured"}")
        println("GCC build script: ${project.findProperty("build.script.gcc")}")
        println("Release build script: ${project.findProperty("build.script.release")}")
        println("testing enabled: ${project.findProperty("test.enabled")}")
        println("test script: ${project.findProperty("test.script") ?: "not configured"}")
        println("test command: ${project.findProperty("test.command") ?: "not configured"}")
    }
}

tasks.named("clean") {
    description = "cleanup build artifacts"
    group = "ci"
    dependsOn("cleanWorkspace")
}