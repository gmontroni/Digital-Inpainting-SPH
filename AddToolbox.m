getd = @(p)path(p,path); % scilab users must *not* execute this

%% add the toolboxes to the path.

getd('packages/');
getd('geodesic/');
getd('geodesic/matlab');
getd('geodesic/src');
getd('toolbox_G/');
getd('toolbox_graph/');
getd('toolbox_graph/mex/');
getd('toolbox_graph/toolbox/');
getd('data/');
getd('save/');

% mesh = mesh3.TriMesh.readPLY('bunny.ply');
% h = mesh.render(axes);
% h.ButtonDownFcn = @buttonDownCallback; % h é o "handle" do "patch"

