
        options(warn=-1)
        package_repo <- 'stephens999/ashr'
        package <- basename(package_repo)
        if (require(package, character.only=TRUE, quietly=TRUE)) {
            write(paste(package, packageVersion(package), "AVAILABLE"), file="/tmp/tmp3mc1svr1.txt")
        } else {
            devtools::install_github(package_repo)
            # if it still does not exist, write the package name to output
            if (require(package, character.only=TRUE, quietly=TRUE)) {
                write(paste(package, packageVersion(package), "INSTALLED"), file="/tmp/tmp3mc1svr1.txt")
            } else {
                write(paste(package, "NA", "MISSING"), file="/tmp/tmp3mc1svr1.txt")
                quit("no")
            }
        }
        cur_version <- packageVersion(package)
        