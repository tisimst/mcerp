Changes
=======

This is the changelog for ``mcerp``.
You can find the releases at https://pypi.org/project/mcerp/

v0.12 (Jan 9, 2019)
-------------------

- Python 3 support added (now ``mcerp`` works on Python 2.7 as well as
  Python 3.5 or later with a single codebase)
- Several small issues fixed

v0.11 (Jun 12, 2014)
--------------------

(changelog wasn't filled, not clear what changed)

v0.10 (Dec 10, 2013)
--------------------

(changelog wasn't filled, not clear what changed)

v0.9.2 (Aug 17, 2013)
---------------------

- Fixed bug in 'umath.py' where the 'mcpts' member references should have been
  '_mcpts'.

v0.9 (Jul 11, 2013)
-------------------

- Removed dependencies on the 'ad' package since the core functionality didn't
  really depend on it.

- Updated the umath sub-module to only utilize numpy for its mathematical
  functions (i.e., sin, exp, abs, etc.).

- Added many constructor functions to make creating UncertainVariables easier,
  like 'Normal', 'Uniform', 'Gamma', etc. All examples have been updated
  accordingly to show the use of these constructors. The original constructor,
  like 'uv(ss.norm(loc=2, scale=0.1))', is still functional.

- Updated the 'describe' method to allow the user to specify a printed 
  identifier as an input other than the 'tag' member.

- Added 'covariance_matrix' and 'correlation_matrix' utility functions.

- Renamed display text from 'UF(...)' and 'UV(...)' to just 'uv(...)'.

- Made the sample statistics permanent properties via the @property decorator.

v0.8 (Jun 27, 2013)
-------------------

- First public release.
 
