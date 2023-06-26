# Decorators for fitting class
# See http://127.0.0.1:8888/lab/tree/tests_demos/decorators_wrapt_tests_140921.ipynb
# 14/09/21

import wrapt

# Wrapt example
def with_arguments(myarg1, myarg2):
    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        return wrapped(*args, **kwargs)
    return wrapper

@with_arguments(1, 2)
def function():
    pass
