from skyfield.api import Star, load
from skyfield.data import hipparcos
from skyfield.framelib import itrs
planets = load('de421.bsp')
earth = planets['earth']

with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

barnards_star = Star.from_dataframe(df.loc[87937])



df = df[df['ra_degrees'].notnull()]
df = df[df['magnitude'] <= 1]

print('After filtering, there are {} stars'.format(len(df)))




ts = load.timescale()
t = ts.now()
bright_stars = Star.from_dataframe(df)
t = ts.utc(2000, 1, 1)
astrometric = earth.at(t).observe(bright_stars)
ra, dec, distance = astrometric.radec()
x, y, z = astrometric.frame_xyz(itrs).au

print(x)