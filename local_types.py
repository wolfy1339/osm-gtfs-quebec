# pylint: disable=missing-module-docstring, missing-class-docstring
from typing import Literal, NotRequired, TypedDict, TypeVar
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

GTFSStopTags = TypedDict('GTFSStopTags', {
    "name": str,
    "ref": str,
    "network": Literal["RTC"],
    "network:wikidata": Literal["Q3456768"],
    "operator": Literal["Réseau de transport de la Capitale"],
    "public_transport": Literal["platform"],
    "highway": Literal["bus_stop"],
    "bus": Literal["yes"],
    "description": NotRequired[str],
    "wheelchair": Literal["yes", "no"],
    "local_ref": NotRequired[str]
})
class GTFSStop(TypedDict):
    props: dict[str, str | int]
    tags: GTFSStopTags
    gtfs_props: dict[str, str]
    geom: Geometry

class RouteRelationMembersPropsWay(TypedDict):
    type: Literal["way"]
    role: str
    ref: str
class RouteRelationMembersPropsNode(TypedDict):
    type: Literal["node"]
    role: str
    ref: str
class RouteRelationMembersPropsRelation(TypedDict):
    ref: int
    type: Literal["relation"]
    role: Literal[""]
class RouteRelationMembers(TypedDict):
    props: RouteRelationMembersPropsWay | RouteRelationMembersPropsNode | RouteRelationMembersPropsRelation
class RouteMemberRelation(TypedDict):
    props: RouteRelationMembersPropsRelation
class RouteRelation(TypedDict):
    props: dict[str, str | int]
    tags: dict[str, str]
    members: dict[str, list[RouteRelationMembers]]

RouteMasterRelationTags = TypedDict('RouteMasterRelationTags', {
    "name": str,
    "ref": str,
    "network": Literal["RTC"],
    "network:wikidata": Literal["Q3456768"],
    "operator": Literal["Réseau de transport de la Capitale"],
    "type": Literal["route_master"],
    "route_master": Literal["bus"],
    "public_transport:version": Literal[2],
    "school": NotRequired[Literal["yes"]],
    "gtfs:route_id": str
})
class RouteMasterRelation(TypedDict):
    props: dict[Literal["id"], int]
    tags: RouteMasterRelationTags
    members: list[RouteMemberRelation]

class RouteStopsData(TypedDict):
    first_stop_name: str
    last_stop_name: str
    trip_id: str
    stops: list[tuple[str, str, str]]

class ExistingStops(TypedDict):
    tags: dict[str, str]
    props: dict[str, str]
    geom: Geometry

class ExistingRouteMasterRelation(TypedDict):
    props: dict[Literal["id"], int]
    tags: dict[str, str]
    members: dict[Literal["relations"], list[RouteMemberRelation]]

class ExistingData(TypedDict):
    stops: list[ExistingStops]
    routes: list[RouteRelation]
    route_masters: list[ExistingRouteMasterRelation]
