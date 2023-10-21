import requests
import pytest
import sys

def main():
    print("Hello world")
    print("Hello world", file=sys.stderr)

if __name__=="__main__":
    main()
