# build a custom exception to handle bad responses from BOLD
class BadResponseError(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class APIOverload(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class ProxyNotWorking(Exception):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
