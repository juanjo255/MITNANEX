"""Allow ``python -m mitnanex ...``."""

import sys

from .cli import main

sys.exit(main())
