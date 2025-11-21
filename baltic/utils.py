import datetime as dt
import re
from itertools import takewhile
from typing import TypeVar, Callable


def decimalDate(date: str, fmt: str = "%Y-%m-%d", variable: bool = False):
    """
    Converts calendar dates in specified format to decimal date.

    A decimal date represents the fraction of the year that has passed by the given date.

    Parameters:
    date (str): The date to be converted.
    fmt (str): The format of the input date string. Default is "%Y-%m-%d".
    variable (bool): If True, allows for variable date precision. Default is False.

    Returns:
    float: The decimal representation of the date.

    Notes:
    - If `variable` is True, the function adjusts the format to match the available precision in the date.
    - For example, a date like "2023" will be interpreted as 2023 Jan 01, while "2023-05" will be interpreted as 2023 May 01.

    Examples:
    >>> decimalDate("2023-05-23")
    2023.3890410958904
    >>> decimalDate("2023", fmt="%Y", variable=True)
    2023.0

    Docstring generated with ChatGPT 4o.
    """
    if not fmt:
        return date

    if variable:  ## if date is variable - extract what is available
        parts = re.split(r"((?<!%)%[YybBmdjWU])", fmt)

        def has_year(s):
            return "%y" in s or "%Y" in s

        formats = [
            # left-to-right Y-m-d -> Y-m -> Y
            *takewhile(has_year, ["".join(parts[:i]) for i in range(len(parts), 0, -1)]),
            # right-to-left d-m-Y -> m-Y -> Y
            *takewhile(has_year, ["".join(parts[i:]) for i in range(len(parts))]),
        ]
        for _fmt in formats:
            try:
                adatetime = dt.datetime.strptime(date, _fmt)
                break
            except ValueError:
                pass
        else:
            adatetime = dt.datetime.strptime(date, fmt)
    else:
        adatetime = dt.datetime.strptime(date, fmt)  ## convert to datetime object

    year = adatetime.year  ## get year
    boy = dt.datetime(year, 1, 1)  ## get beginning of the year
    eoy = dt.datetime(year + 1, 1, 1)  ## get beginning of next year
    return year + ((adatetime - boy).total_seconds() / ((eoy - boy).total_seconds()))  ## return fractional year


def calendarDate(timepoint: float, fmt: str = "%Y-%m-%d"):
    """
    Converts decimal dates to a specified calendar date format.

    A decimal date represents the fraction of the year that has passed by the given timepoint.
    This function converts it back to a calendar date in the given format.

    Parameters:
    timepoint (float): The decimal representation of the date.
    fmt (str): The desired format of the output date string. Default is '%Y-%m-%d'.

    Returns:
    str: The date in the specified calendar format.

    Examples:
    >>> calendarDate(2023.3923497267758)
    '2023-05-24'
    >>> calendarDate(2023.0, fmt='%Y')
    '2023'

    Docstring generated with ChatGPT 4o.
    """
    year = int(timepoint)
    rem = timepoint - year

    base = dt.datetime(year, 1, 1)
    result = base + dt.timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)

    return dt.datetime.strftime(result, fmt)


def convertDate(date_string: str, start: str, end: str):
    """
    Converts calendar dates between given formats.

    Parameters:
    x (str): The date string to be converted.
    start (str): The format of the input date string.
    end (str): The desired format of the output date string.

    Returns:
    str: The date converted to the new format.

    Examples:
    >>> convertDate('23-05-2023', '%d-%m-%Y', '%Y/%m/%d')
    '2023/05/23'
    >>> convertDate('2023/05/23', '%Y/%m/%d', '%B %d, %Y')
    'May 23, 2023'

    Docstring generated with ChatGPT 4o.
    """
    return dt.datetime.strftime(dt.datetime.strptime(date_string, start), end)
    try:
        date_obj = dt.datetime.strptime(date_string, start)
        return dt.datetime.strftime(date_obj, end)
    except ValueError as e:
        raise ValueError('Error converting date "%s" from format "%s" to "%s": "%s"' % (date_string, start, end, e))


S = TypeVar("S")
T = TypeVar("T")


def initialized_property(func: Callable[[S], T]):
    name = func.__name__

    def getter(self: S) -> T:
        value = getattr(self, f"_{name}", None)
        if value is None:
            raise AttributeError(f"{name} is not initialized")
        return value

    def setter(self: S, value: T):
        setattr(self, f"_{name}", value)

    return property(getter, setter)


def always_true(*args):
    return True
