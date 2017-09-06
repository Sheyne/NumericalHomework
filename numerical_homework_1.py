import math

R = 6367444.50 # earth radius in meters
s = 86164.09 # sidereal day in seconds
h = 20200000 # satelite altitude in m
p = s/2 # period of satelite orbit in seconds
c = 3e8

should_wrap_ranges = True

def wrap_range(n, low, high):
	if not should_wrap_ranges:
		return n
	return (n - low) % (high - low) + low

def test_wrap_range():
	assert wrap_range(3, -2, 2) == -1
	assert wrap_range(-2, -2, 2) == -2
	assert wrap_range(2, -2, 2) == -2
	assert wrap_range(4, -2, 2) == 0
	assert wrap_range(1, -2, 2) == 1

def dms_to_radians(degrees, minutes, seconds, sgn=1):
	"convert degrees minutes seconds to radians"

	# Pretty straightforward, sums degrees, minutes / 60, and seconds / 3600 to get total
	# degrees then multiplies by pi / 180 to convert to radians. Finally multiply in the sign
	# once we've determined the magnitude of rads
	
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
	"convert radians to degrees minutes seconds"

	# get the decimal degrees from the radians
	deg_raw = (rad * 180 / math.pi)

	# note sign and enforce positivity
	sgn = 1 if deg_raw >= 0 else -1
	deg_raw *= sgn

	# floor decimal degrees to get degrees
	deg = int(deg_raw)
	# floor decimal minutes to get minutes
	min_raw = (deg_raw - deg) * 60
	min_ = int(min_raw)
	# round seconds to get seconds
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
	"""
	Get the position from lat, lon, alt to cartesian, assuming the earth isn't spinning
	"""

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

def position_of_meridian(t):
	return t * 360 / s 

def lat_lon_to_cart(t, lat, lon, alt):
	"""
	Add the earth's rotation to the longitude then convert to cartesian. 
	"""
	lons_to_add = position_of_meridian(t)

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


def above_horizon(earth_point, sat_point):
	"""
	Determine if the satellite is above the horizon. Effectively we want to 
	get the scalar projection of the satellite vector onto the earth position vector. 

	scalar projection of a onto b is 

	p = (|a| / |b|) cos(th)

	a.b = |a||b| cos(th)
	cos(th) = a.b / (|a||b|)

	p = (|a| / |b|) a.b / (|a||b|)
	p = a.b / |b|^2

	is this projection scales the earth vector by more than 1, the satellite is 
	above the horizon, otherwise the satellite is below the horizon
	"""

	return (earth_point @ sat_point) / (sat_point ** 2).sum() >= 1

def timeless_cart_to_lat_lon(v):
	x, y, z = v
	r = (v ** 2).sum() ** 0.5
	alt = r - R

	lon = math.atan2(y, x)
	lat = math.asin(z / r)

	return radians_to_dms(lat), radians_to_dms(lon), alt


def test_timeless_cart_to_lat_lon():
	import numpy
	import pytest

	def test(lat_in, lon_in, alt_in):
		cart = numpy.array(timeless_lat_lon_to_cart(lat_in, lon_in, alt_in))
		print("cart:", cart)
		lat, lon, alt = timeless_cart_to_lat_lon(cart)
		assert lat == lat_in
		assert lon == lon_in
		assert pytest.approx(alt) == alt_in

	test((45, 0, 0, 1), (45, 0, 0, 1), 200)
	test((0, 0, 0, 1), (45, 0, 0, 1), 200)
	test((0, 0, 0, 1), (0, 0, 0, 1), 200)
	test((0, 0, 0, 1), (90, 0, 0, 1), 200)
	test((90, 0, 0, 1), (90, 0, 0, 1), 200)


def cart_to_lat_lon(t, v):
	lat, lon, alt = timeless_cart_to_lat_lon(v)

	return lat, lon - position_of_meridian(t), alt
