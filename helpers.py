import numpy as np
import time
import random
from astropy.time import Time

def find_oc_outliers(ephemeris_obj, y_limit, return_type='str'):
    # Used in O-C plotting to calculate scatter points
    DAYS_TO_SECONDS = 86400
    lin_model = ephemeris_obj.get_model_ephemeris('linear')
    y = (ephemeris_obj._subtract_plotting_parameters(ephemeris_obj.timing_data.mid_times, lin_model['conjunction_time'], lin_model['period'], ephemeris_obj.timing_data.epochs)) * DAYS_TO_SECONDS
    outliers = []
    outliers_str = []
    i = 0
    while i in range(len(y)):
        # if epochs[i] > 3000 and epochs[i] < 4000:
        if y[i] < y_limit:
            outliers.append([i, ephemeris_obj.timing_data.epochs[i], y[i], ephemeris_obj.timing_data.mid_times[i]])
            outliers_str.append(f"index: {i}, epoch: {ephemeris_obj.timing_data.epochs[i]}, y value: {y[i]}, mid-time: {ephemeris_obj.timing_data.mid_times[i]}")
        i+=1
    if return_type == 'str':
        return outliers_str
    elif return_type == 'numerical':
        return outliers
    
def convert_to_tdb(skycoord_obj, location_obj, time_obj):
    """
    time_obj:
        Astropy time object, must have time format of JD, can be any time scale.
    """
    # Set location for time object
    time_obj.location = location_obj
    # Calc barycentric light travel time
    ltt_bary = time_obj.light_travel_time(skycoord_obj)
    # Add light travel time and set scale to TDB
    times_bjd = time_obj.tdb + ltt_bary
    return times_bjd

def get_epochs_for_new_data(T, T0, P, tra_or_occ):
    N = (T-T0)/P
    if tra_or_occ == "occ":
        N = np.floor(N)
    return int(N)
    
def str_time_prop(start, end, time_format, prop):
    """Get a time at a proportion of a range of two formatted times.

    start and end should be strings specifying times formatted in the
    given format (strftime-style), giving an interval [start, end].
    prop specifies how a proportion of the interval to be taken after
    start.  The returned time will be in the specified format.
    """

    stime = time.mktime(time.strptime(start, time_format))
    etime = time.mktime(time.strptime(end, time_format))

    ptime = stime + prop * (etime - stime)

    return time.strftime(time_format, time.localtime(ptime))


def random_date(start, end, prop):
    return str_time_prop(start, end, '%Y-%m-%d %H:%M', prop)

def get_n_random_observations(n, tra_or_occ, start_time, end_time, astroplan_target, skycoord_obj, earthlocation_obj, conjunction_time, orb_period):
    random.seed(1) # Set seed so we get the same random dates (this ensures reproducibility)
    # Get n random dates
    random_dates = []
    for i in range(n):
        random_dates.append(random_date(start_time, end_time, random.random()))
    # Get next occultation for each date
    new_mid_times = []
    for date in random_dates:
        if tra_or_occ == "tra":
            new_mid_times.append(astroplan_target.next_primary_eclipse_time(Time(date), n_eclipses=1)[0])
        elif tra_or_occ == "occ":
            new_mid_times.append(astroplan_target.next_secondary_eclipse_time(Time(date), n_eclipses=1)[0])
    # Convert list to one Time object and convert to JD format
    new_mid_times = Time(new_mid_times)
    new_mid_times.format = "jd"
    # Convert each occultation into BJD TDB
    new_mid_times_bjd = convert_to_tdb(skycoord_obj, earthlocation_obj, new_mid_times)
    # Get epoch for each time
    new_epochs = []
    for time in new_mid_times_bjd:
        new_epochs.append(get_epochs_for_new_data(np.array([time.value]), conjunction_time, orb_period, tra_or_occ))
    new_errs = np.full(n, 0.0001)
    new_tra_or_occ = np.full(n, "occ")
    return_dict = {
        "new_epochs": new_epochs,
        "new_mid_times": new_mid_times, 
        "new_mid_time_errs": new_errs, 
        "new_tra_or_occ": new_tra_or_occ
    }
    return return_dict