import os
import sys

# liftover/ scripts import siblings as top-level modules (as when run as scripts)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'liftover'))
