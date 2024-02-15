
# MotifFindR

<!-- badges: start -->
<!-- badges: end -->

`MotifFindR` is an R project used to identify Transcription Factor Binding Sites through motif analysis and was created exclusively for an internal renal cancer research project, 2023-2024.

## Installation

- **Install Git:** If you haven't already, download and install Git from [git-scm.com](https://git-scm.com/).
- **Open a Terminal or Command Prompt:** Access your command line interface. On Windows, you can use Command Prompt or PowerShell. On macOS or Linux, open the Terminal.
- **Navigate to the Directory:** Use the `cd` command to navigate to the directory where you want to clone the repository.
- **Clone the Repository:** Use the git clone command followed by the URL of the repository.  See below:

``` 
git clone https://github.com/timyers/MotifFindR
```
- `MotifFindR` uses the `R` package [`renv`](https://rstudio.github.io/renv/index.html) to create a reproducible environment.  It records the exact package versions `MotifFindR` depends on in a lockfile, *renv.lock*, and ensures those exact versions get installed.

- After installation, when you open the project `renv` should automatically bootstrap itself, downloading and installing the appropriate version of `renv`.  It should also ask if you want to download and install all the packages `MotifFindR` needs by running `renv::restore()`.

- **NOTE:** The `renv` package tracks but does not help with the version of R used (the version of `R` used can be found in the `renv.lock` file).  `renv` can't help with this because it runs inside of `R`.  However, there are other tools that can help.  But tools like [The R Installation Manager or `rig`](https://github.com/r-lib/rig) can help with switching between multiple versions of `R` on one computer.

## Example

This is a basic example.

``` r
## basic example code
```

