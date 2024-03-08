from typing import Literal, Optional, TypedDict, TypeVar
from osgeo.ogr import Geometry

K = TypeVar('K')
class GTFSDataKey[K](TypedDict):
    field_names: list[K]
    data: list[dict[K, str]]

class GTFSData(TypedDict):
    agency: GTFSDataKey[str]
    calendar_dates: GTFSDataKey[str]
    feed_info: GTFSDataKey[str]
    Horaire_Boucle_Partage: GTFSDataKey[str]
    routes: GTFSDataKey[str]
    shapes: GTFSDataKey[str]
    stop_times: GTFSDataKey[str]
    stops: GTFSDataKey[str]
    transfers: GTFSDataKey[str]
    trips: GTFSDataKey[Literal['route_id', 'service_id', 'shape_id', 'trip_id', 'trip_headsign', 'trip_short_name', 'direction_id', 'block_id','wheelchair_accessible']]

class GTFSStopTags(TypedDict):
    bus: str
    highway: str
    name: str
    description: str
    public_transport: str
    ref: str
    network: str
    operator: str
    wheelchair: str
    local_ref: Optional[str]

class GTFSStop(TypedDict):
    props: dict[str, str]
    tags: GTFSStopTags
    gtfs_props: dict[str, str]
    geom: Geometry
