#' 
#' @title Universal Thermal Climate Index in climate4R
#' @description Computation of the Universal Thermal Climate Index (UTCI) by ECA&D (European Climate 
#' Assessment & Dataset) directly from climate4R objects.
#'
#' @param tmax A climate4R dataset of daily maximum air temperature.
#' @param tmin A climate4R dataset of daily minimum air temperature.
#' @param hurs A climate4R dataset of daily relative humidity.
#' @param wind A climate4R dataset of daily wind speed.
#' @param radiation A climate4R dataset of daily radiation.
#' @param mask Optional. A binary climate4R grid (0 and 1) with lat x lon dimensions to mask the output. 
#' UTCI will be computed at gridboxes with mask value 1, while the rest will be set to NA.
#' If NULL (default), UTCI is calculated at all gridboxes (sea values may be inaccurate).
#' @param parallel Logical. Enable parallel processing (default = FALSE).
#' @param max.ncores Integer. Maximum number of cores to use (default = 16).
#' @param ncores Integer. Specific number of cores to use.
#'
#' @return A climate4R object with the computed index.
#' @details Air temperature datasets are internally converted to degC, while relative humidity is converted to %,
#' wind speed to m s-1 and radiation to MJ m-2. 
#' The Fortran backend \code{calc_UTCI_vector} of the \pkg{UTCIr} package is internally called.
#'
#' @references
#' Bröde, P. et al. (2012). Deriving the operational procedure for the Universal Thermal
#' Climate Index (UTCI). \emph{International Journal of Biometeorology}, 56, 481-494.
#' \doi{10.1007/s00484-011-0454-1}
#'
#' @import transformeR
#' @importFrom magrittr %>% %<>%
#' @importFrom convertR udConvertGrid
#' @importFrom udunits2 ud.are.convertible ud.convert
#' @importFrom parallel stopCluster
#' 
#' @examples \dontrun{
#' require(climate4R.UTCI)
#'
#' data("ERA5_day_t2mx", package = "climate4R.UTCI")
#' data("ERA5_day_t2mn", package = "climate4R.UTCI")
#' data("ERA5_day_hurs", package = "climate4R.UTCI")
#' data("ERA5_day_sfcwind", package = "climate4R.UTCI")
#' data("ERA5_day_ssrd", package = "climate4R.UTCI")
#'
#' utci <- utciGrid(tmax = ERA5_day_t2mx, tmin = ERA5_day_t2mn,
#'                  hurs = ERA5_day_hurs, wind = ERA5_day_sfcwind,
#'                  radiation = ERA5_day_ssrd)
#' }
#'
#' @author climate4R adaptation by C. Rodriguez-Rumayor.
#' Original \pkg{UTCIr} by ECA-D (\url{https://github.com/ECA-D/UTCIr}).
#' @export

utciGrid <- function(tmax,
                     tmin,
                     hurs,
                     wind,
                     radiation,
                     mask = NULL,
                     parallel = FALSE,
                     max.ncores = 16,
                     ncores = NULL) {
    
    # Basic checks
    missing.vars <- c(
        if (is.null(tmax)) "tmax",
        if (is.null(tmin)) "tmin",
        if (is.null(hurs)) "hurs",
        if (is.null(wind)) "wind",
        if (is.null(radiation)) "radiation")
    if (length(missing.vars) > 0) {
        stop("Some input variables are missing. ", paste(missing.vars, collapse = ", "), " required for UTCI calculation.")
    }

    stopifnot(all(sapply(list(tmax, tmin, hurs, wind, radiation), function(x) isGrid(x))))
    if (any(sapply(list(tmax, tmin, hurs, wind, radiation), isMultigrid))) {
        stop("Multigrids are not an allowed input")
    }

    grid_types <- sapply(list(tmax, tmin, hurs, wind, radiation), typeofGrid)
    if (length(unique(grid_types)) > 1) {
        stop("All input variables must be of the same type (either grid or station).")
    }

    if(is.null(mask)) message("NOTE: No mask provided. UTCI values over sea may be inaccurate.")

    # Convert inputs to required units (if needed)
    for (temp in c("tmax", "tmin")) {
        obj <- get(temp)
        temp.u <- getGridUnits(obj)
        if (tolower(temp.u) %in% c("degrees celsius", "degree celsius")) {
            attr(obj$Variable, "units") <- "degC"
            temp.u <- "degC"
        }
        if (tolower(temp.u) %in% c("degrees fahrenheit", "degree fahrenheit")) {
            attr(obj$Variable, "units") <- "degF"
            temp.u <- "degF"
        }
        if (ud.are.convertible(temp.u, "degC")) {
            if (ud.convert(1, temp.u, "degC") != 1) {
                message("[", Sys.time(), "] Converting ", temp, " units ...")
                obj %<>% udConvertGrid(new.units = "degC")
            }
        } else {
            stop("Non compliant ", temp, " units (", temp.u, " is not convertible to degC)")
        }
        assign(temp, obj, envir = environment())
    }
    hurs.u <- getGridUnits(hurs)
    if (tolower(hurs.u) %in% c("percentage")) {
        attr(hurs$Variable, "units") <- "%"
        hurs.u <- "%"
    }
    if (ud.are.convertible(hurs.u, "%")) {
        if (ud.convert(1, hurs.u, "%") != 1) { 
            message("[", Sys.time(), "] Converting relative humidity units ...")
            hurs %<>% udConvertGrid(new.units = "%") 
        }  
    } else {
        stop("Non compliant hurs units (", hurs.u, " is not convertible to %)")
    }
    wind.u <- getGridUnits(wind)
    if (ud.are.convertible(wind.u, "m s-1")) {
        if (ud.convert(1, wind.u, "m s-1") != 1) { 
            message("[", Sys.time(), "] Converting wind speed units ...")
            wind %<>% udConvertGrid(new.units = "m s-1") 
        }  
    } else {
        stop("Non compliant wind units (", wind.u, " is not convertible to m s-1)")
    }
    # Allow conversion of radiation from energy flux to power flux (e.g. J m-2 to W m-2) 
    rad.u <- getGridUnits(radiation)
    if (ud.are.convertible(rad.u, "MJ m-2")) {
        if (ud.convert(1, rad.u, "MJ m-2") != 1) {
            message("[", Sys.time(), "] Converting radiation units ...")
            radiation %<>% udConvertGrid(new.units = "MJ m-2")
        } 
    } else if (ud.are.convertible(rad.u, "W m-2")) {
        if (ud.convert(1, rad.u, "W m-2") != 1) {
            radiation %<>% udConvertGrid(new.units = "W m-2")
        }
        message("[", Sys.time(), "] Converting radiation units ...")
        time_step <- difftime(as.POSIXct(getRefDates(radiation)[2], tz = "UTC"),
                              as.POSIXct(getRefDates(radiation)[1], tz = "UTC"),
                              units = "secs") %>% as.numeric()
        radiation$Data <- radiation$Data * time_step / 1e6
        attr(radiation$Variable, "units") <- "MJ m-2"
    } else {
        stop("Non compliant radiation units (", rad.u, " is not convertible to 'W m-2' or 'MJ m-2')")
    } 

    # Sanity check on hurs values
    if (any(hurs$Data < 0 | hurs$Data > 100, na.rm = TRUE)) {
        stop("Some relative humidity values are outside the expected [0, 100] range")
    }

    # Replace possible Inf/-Inf values with NA
    for (x in c("tmax", "tmin", "hurs", "wind", "radiation")) {
        obj <- get(x)
        obj$Data[is.infinite(obj$Data)] <- NA
        assign(x, obj, envir = environment())
    }

    # Ensure member dimension
    tmax %<>% redim(member = TRUE)
    tmin %<>% redim(member = TRUE)
    hurs %<>% redim(member = TRUE)
    wind %<>% redim(member = TRUE)
    radiation %<>% redim(member = TRUE)

    # Consistency checks
    suppressMessages(checkDim(tmax, tmin, hurs, wind, radiation, dimensions = c("time", "lat", "lon")))

    station <- typeofGrid(tmax) == "station"
    n.mem <- getShape(tmax, "member")

    # Parallel processing setup
    if (n.mem > 1) {
        parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
        apply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "lapply")
        if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
    } else {
        if (isTRUE(parallel)) message("NOTE: Parallel processing was skipped (unable to parallelize one single member)")
        apply_fun <- lapply
    }

    # Coordinates and dates
    dates <- as.Date(getRefDates(tmax))
    dates.jd <- as.integer(format(dates, "%j"))
    coords <- getCoordinates(tmax)

    if (station) {
        lats <- coords$y
        n.loc <- length(lats)
    } else {
        # For grids, build flat vectors of lat per cell
        lats <- rep(coords$y, times = length(coords$x))
        n.loc <- length(lats)
    }

    # Build ind of locations/gridboxes to compute UTCI on
    if (!is.null(mask)) {
        # Drop residual dimensions and check binary nature
        mask$Data <- drop(mask$Data)
        vals <- unique(as.vector(mask$Data))
        if (!all(vals[!is.na(vals)] %in% c(0, 1))) stop("Provided mask is not a binary grid.")
    
        if (station) {
            # For each station, find the nearest mask gridbox
            mask.stations <- interpGrid(mask,
                                        new.coordinates = getCoordinates(tmax),
                                        method = "nearest")
            ind <- which(as.numeric(mask.stations$Data) > 0)
        } else {
            # For grids, check consistency with input variables
            attr(mask$Data, "dimensions") <- c("lat", "lon")
            tryCatch(
                suppressMessages(checkDim(tmax, mask, dimensions = c("lat", "lon"))),
                error = function(e) stop("Provided mask and input variables have a different spatial resolution.", call. = FALSE)
            )
            ind <- which(as.numeric(mask$Data) > 0)
        }
    } else {
        ind <- seq_len(n.loc)
    }

    # Load Fortran backend
    dyn.load(system.file("libs/UTCIr.so", package = "UTCIr"))

    message("[", Sys.time(), "] Calculating UTCI ...")

    # Compute UTCI member by member
    out.list <- apply_fun(seq_len(n.mem), function(m) {

        tmax_m <- subsetGrid(tmax, members = m, drop = TRUE)[["Data"]]; if (!station) tmax_m <- array3Dto2Dmat(tmax_m)
        tmin_m <- subsetGrid(tmin, members = m, drop = TRUE)[["Data"]]; if (!station) tmin_m <- array3Dto2Dmat(tmin_m)
        hurs_m <- subsetGrid(hurs, members = m, drop = TRUE)[["Data"]]; if (!station) hurs_m <- array3Dto2Dmat(hurs_m)
        wind_m <- subsetGrid(wind, members = m, drop = TRUE)[["Data"]]; if (!station) wind_m <- array3Dto2Dmat(wind_m)
        rad_m <- subsetGrid(radiation, members = m, drop = TRUE)[["Data"]]; if (!station) rad_m <- array3Dto2Dmat(rad_m)

        n.time <- nrow(tmax_m)
        utci.mat <- matrix(NA, nrow = n.time, ncol = n.loc)

        for (i in ind) {

            tmax.col <- tmax_m[, i]; tmax.col[is.na(tmax.col)] <- -999.9
            tmin.col <- tmin_m[, i]; tmin.col[is.na(tmin.col)] <- -999.9
            hurs.col <- hurs_m[, i]; hurs.col[is.na(hurs.col)] <- -999.9
            wind.col <- wind_m[, i]; wind.col[is.na(wind.col)] <- -999.9
            rad.col  <- rad_m[, i];  rad.col[is.na(rad.col)]   <- -999.9

            utci.vec <- as.double(rep(-999.9, n.time))
            result <- .Fortran("calc_UTCI_vector",
                               dates.jd,
                               as.double(lats[i]),
                               as.double(tmax.col),
                               as.double(tmin.col),
                               as.double(hurs.col),
                               as.double(rad.col),
                               as.double(wind.col),
                               as.integer(n.time),
                               utci.vec,
                               PACKAGE = "UTCIr")
            utci.vec <- result[[9]]
            utci.vec[utci.vec < -999.9] <- NA
            utci.mat[, i] <- utci.vec
        }

        # Build output grid
        out.grid <- subsetGrid(tmax, members = m, drop = FALSE)     
        if (station) {
            out.grid[["Data"]] <- utci.mat
            attr(out.grid[["Data"]], "dimensions") <- c("time", "loc")
        } else {
            xy <- getCoordinates(tmax)
            out.grid[["Data"]] <- mat2Dto3Darray(utci.mat, x = xy$x, y = xy$y)
        }

        # Update variable metadata
        out.grid[["Variable"]] <- list(varName = "UTCI", level = NULL)
        attr(out.grid[["Variable"]], "units") <- "degC"
        attr(out.grid[["Variable"]], "longname") <- "Universal Thermal Climate Index"

        return(out.grid)
    })

    # Combine members
    out <- if (length(out.list) == 1) {
        out.list[[1]]
    } else {
        suppressWarnings(do.call(bindGrid, c(out.list, list(dimension = "member"))))
    }

    # Final redim for station data
    if (station) out %<>% redim(drop = FALSE, loc = TRUE, member = FALSE)

    message("[", Sys.time(), "] Done.")
    invisible(out)

}
