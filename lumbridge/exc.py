class LumbridgeException(BaseException):
    def __init__(self, message: str = ""):
        super().__init__(message)
        self.message = message

    def __getitem__(self):
        return self.message


class WrongInputPath(LumbridgeException):
    def __init__(self, path: str | None = None):
        message = f"The provided path is not a folder. Problem at: {path}"
        super().__init__(message)


class WrongFile(LumbridgeException):
    def __init__(self, path: str, expected_ending: str):
        message = f"The provided file at {path} has a wrong file ending. Expected {expected_ending}"
        super().__init__(message)


class MissingFile(LumbridgeException):
    def __init__(self, missing_files: dict[str, set[str]]):
        message = f"The input has missing files: {missing_files}"
        super().__init__(message)
