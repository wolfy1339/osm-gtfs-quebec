# pylint: disable=missing-module-docstring, missing-class-docstring
from typing import Literal, NotRequired, TypedDict, TypeVar, Generic
from osgeo.ogr import Geometry

K = TypeVar('K')
class GTFSDataKey(TypedDict, Generic[K]):
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
    "network": Literal["RTC"] | Literal["STLévis"],
    "network:wikidata": Literal["Q3456768"] | Literal["Q3488027"],
    "operator": Literal["Réseau de transport de la Capitale"] | Literal["Société de transport de Lévis"],
    "public_transport": Literal["platform"],
    "highway": Literal["bus_stop"],
    "bus": Literal["yes"],
    "description": NotRequired[str],
    "wheelchair": Literal["yes", "no"],
    "local_ref": NotRequired[str],
    "gtfs:stop_id": str
})
class GTFSStop(TypedDict):
    props: dict[str, str | int]
    tags: GTFSStopTags
    gtfs_props: dict[str, str]
    geom: Geometry

Type = TypeVar('Type', bound=Literal["way", "node", "relation"])
Role = TypeVar('Role', bound=str)
class RelationMembersProps(TypedDict, Generic[Type, Role]):
    type: Type
    role: Role
    ref: str
class RouteRelationMemberProps(TypedDict, Generic[Type, Role]):
    props: RelationMembersProps[Type, Role]

RelationMemberNodes = RouteRelationMemberProps[Literal["node"], str]
RelationMemberWays = RouteRelationMemberProps[Literal["way"], str]
RouteMemberRelation = RouteRelationMemberProps[Literal["relation"], Literal[""]]

class RouteRelationMembers(TypedDict):
    nodes: list[RelationMemberNodes]
    ways: NotRequired[list[RelationMemberWays]]

class RouteRelation(TypedDict):
    props: dict[str, str | int]
    tags: dict[str, str]
    members: RouteRelationMembers

RouteMasterRelationTags = TypedDict('RouteMasterRelationTags', {
    "name": str,
    "ref": str,
    "network": Literal["RTC"] | Literal["STLévis"],
    "network:wikidata": Literal["Q3456768"] | Literal["Q3488027"],
    "operator": Literal["Réseau de transport de la Capitale"] | Literal["Société de transport de Lévis"],
    "type": Literal["route_master"],
    "route_master": Literal["bus"],
    "public_transport:version": Literal[2],
    "school": NotRequired[Literal["yes"]],
    "gtfs:route_id": str,
    "colour": NotRequired[str]
})
class RouteMasterRelation(TypedDict):
    props: dict[Literal["id"], int]
    tags: RouteMasterRelationTags
    members: list[RouteMemberRelation]

class RouteStopsData(TypedDict):
    first_stop_name: str
    last_stop_name: str
    trip_id: str
    stops: list[tuple[str | int, str, str]]

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

class StopsData(TypedDict):
    stop_id: str
    stop_name: str
    stop_sequence: int
    pickup_type: str
    drop_off_type: str

class TripsData(TypedDict):
    stops: list[StopsData]
    stop_count: int
    direction_id: str
    trip_headsign: str
