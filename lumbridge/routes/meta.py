from fastapi import APIRouter


meta_router = APIRouter(prefix="/api/meta", tags=["meta"])


@meta_router.get("/info")
async def get_version():
    return {"version": "0.0.1"}


@meta_router.get("/server_ping")
async def get_server_ping():
    return {"server_ping": "pong"}


@meta_router.get("/db_ping")
async def get_db_ping():
    return {"db_ping": "pong"}