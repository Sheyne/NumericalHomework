import math

R = 6367444.50 # earth radius in meters
s = 86164.09 # sidereal day in seconds
h = 20200000 # satelite altitude in m
p = s/2 # period of satelite orbit in seconds

def wrap_range(n, low, high):
	return (n - low) % (high - low) + low

def test_wrap_range():
	assert wrap_range(3, -2, 2) == -1
	assert wrap_range(-2, -2, 2) == -2
	assert wrap_range(2, -2, 2) == -2
	assert wrap_range(4, -2, 2) == 0
	assert wrap_range(1, -2, 2) == 1

def dms_to_radians(degrees, minutes, seconds, sgn=1):
	rads = ((degrees + minutes / 60 + seconds / 60 / 60) * math.pi * sgn / 180)
	return wrap_range(rads, -math.pi, math.pi)


def test_dms_to_radians():
	import pytest
	assert dms_to_radians(0, 0, 0) == pytest.approx(0)
	assert dms_to_radians(360, 0, 0) == pytest.approx(0)
	assert dms_to_radians(365, 0, 0) == pytest.approx(dms_to_radians(5, 0, 0))
	assert dms_to_radians(120, 35, 27) == pytest.approx(2.104707089)
	assert dms_to_radians(5, 16, 45) == pytest.approx(0.09213884)
	assert dms_to_radians(359, 0, 0) == pytest.approx(dms_to_radians(-1, 0, 0))
	assert dms_to_radians(359, 0, 0) < 0
	assert dms_to_radians(160, 0, 0, -1) == pytest.approx(dms_to_radians(200, 0, 0, 1))

def radians_to_dms(rad):
	deg_raw = wrap_range((rad * 180 / math.pi), -180, 180)
	sgn = 1 if deg_raw >= 0 else -1
	deg_raw *= sgn
	deg = int(deg_raw)
	min_raw = (deg_raw - deg) * 60
	min_ = int(min_raw)
	sec = int(round((min_raw - min_) * 60))
	# because of rounding, sec can end up at 60, which causes up to carry
	if sec >= 60:
		sec -= 60
		min_ += 1
	if min_ >= 60:
		min_ -= 60
		deg += 1

	return (deg, min_, sec, sgn)

def test_radians_to_dms():
	assert radians_to_dms(dms_to_radians(120, 35, 27)) == (120, 35, 27, 1)
	assert radians_to_dms(dms_to_radians(367, 35, 27)) == (7, 35, 27, 1)
	assert radians_to_dms(dms_to_radians(7, 0, 0)) == (7, 0, 0, 1)
	assert radians_to_dms(dms_to_radians(7, 1, 1)) == (7, 1, 1, 1)
	assert radians_to_dms(dms_to_radians(7, 1, 2)) == (7, 1, 2, 1)
	assert radians_to_dms(dms_to_radians(183, 0, 0)) == (177, 0, 0, -1)

def timeless_lat_lon_to_cart(lat, lon, alt):
	lat = dms_to_radians(*lat)
	lon = dms_to_radians(*lon)

	r = R + alt
	return (r * math.cos(lon) * math.cos(lat), r * math.sin(lon) * math.cos(lat), r * math.sin(lat))

def test_timeless_lat_lon_to_cart():
	import pytest
	assert timeless_lat_lon_to_cart((0,  0, 0,  1), (  0, 0, 0,  1), 0) == pytest.approx((R, 0, 0) , abs=1e-7)
	assert timeless_lat_lon_to_cart((0,  0, 0,  1), (180, 0, 0,  1), 0) == pytest.approx((-R, 0, 0), abs=1e-7)
	assert timeless_lat_lon_to_cart((0,  0, 0,  1), ( 90, 0, 0,  1), 0) == pytest.approx((0, R, 0) , abs=1e-7)
	assert timeless_lat_lon_to_cart((0,  0, 0,  1), ( 90, 0, 0, -1), 0) == pytest.approx((0, -R, 0), abs=1e-7)
	assert timeless_lat_lon_to_cart((90, 0, 0,  1), (  0, 0, 0,  1), 0) == pytest.approx((0, 0, R) , abs=1e-7)
	assert timeless_lat_lon_to_cart((90, 0, 0, -1), (  0, 0, 0,  1), 0) == pytest.approx((0, 0, -R), abs=1e-7)
	assert timeless_lat_lon_to_cart((90, 0, 0, -1), ( 40, 0, 0,  1), 0) == pytest.approx((0, 0, -R) , abs=1e-7)
	assert timeless_lat_lon_to_cart((45, 0, 0,  1), (  0, 0, 0,  1), 0) == pytest.approx((R * math.sqrt(2)/2, 0, R * math.sqrt(2) / 2) , abs=1e-7)
	assert timeless_lat_lon_to_cart((0, 0, 0, 1), (0, 0, 0, 1), 500) == pytest.approx((R + 500, 0, 0) , abs=1e-7)

def lat_lon_to_cart(t, lat, lon, alt):
	lons_to_add = t * 360 / s 

	d, *a = lon
	lon = (d + lons_to_add, *a)

	return timeless_lat_lon_to_cart(lat, lon, alt)

def test_lat_lon_to_cart():
	import pytest
	assert lat_lon_to_cart(s / 4 + s, (0, 0, 0, 1), (0, 0, 0,  1), 0) == pytest.approx((0, R, 0) , abs=1e-7)
	assert lat_lon_to_cart(s, (0, 0, 0, 1), (0, 0, 0,  1), 0) == pytest.approx((R, 0, 0) , abs=1e-7)
	assert lat_lon_to_cart(s / 2, (0, 0, 0, 1), (0, 0, 0,  1), 0) == pytest.approx((-R, 0, 0) , abs=1e-7)
	assert lat_lon_to_cart(3 * s / 4, (0, 0, 0, 1), (0, 0, 0,  1), 0) == pytest.approx((0, -R, 0) , abs=1e-7)

def sat_pos(t, u, v, alt, th):
	arg = 2 * math.pi * t / p + th
	return (R + alt) * (u*math.cos(arg) + v*math.sin(arg))

def test_sat_pos():
	import numpy
	import pytest

	a = 2.020000000000000000E+07
	u = numpy.array([1,0,0])
	v = numpy.array([0,.57357643,.81915204])
	assert sat_pos(0, u, v, a, 0) == pytest.approx([R + a, 0, 0])
	assert sat_pos(p/2, u, v, a, 0) == pytest.approx([-R - a, 0, 0], abs=1e-8)
	assert sat_pos(p / 4, u, v, a, 0) == pytest.approx(v * (R + a), abs=1e-8)
