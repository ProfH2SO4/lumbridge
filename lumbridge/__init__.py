from fastapi import FastAPI

from lumbridge.routes import meta_router, user_router


__version__ = "0.0.1"
__started_at__ = "2023-02-25T20:00:00"


class Lumbridge(FastAPI):
    def __init__(self):
        super().__init__()



app = FastAPI()


app.include_router(meta_router)
app.include_router(user_router)

