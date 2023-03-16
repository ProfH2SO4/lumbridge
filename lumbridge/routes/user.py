from fastapi import APIRouter
from fastapi.responses import JSONResponse

from lumbridge.enum_ import HTTPStatus

user_router = APIRouter(prefix="/api/users", tags=["users"])



@user_router.put("/{user_id}/create")
async def create_user(user_id: str):
    _ = user_id
    return JSONResponse(status_code=HTTPStatus.CREATED, content={"Status": 1})