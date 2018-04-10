% Math 226B - Homework #2
% Problem 2a
% You are given an 9 × 9 matrix A ? 0 with the following sparsity structure
% Show the graph associated with A.

  % 1,2,3,4,5,6,7,8,9
A =[1,0,0,1,0,0,0,0,1; %1
    0,1,0,1,1,1,1,0,0;  %2
    0,0,1,0,0,0,0,1,0;  %3
    1,1,0,1,0,1,0,0,0;  %4
    0,1,0,0,1,0,0,1,0;  %5
    0,1,0,1,0,1,1,0,1;  %6
    0,1,0,0,0,1,1,1,0;  %7
    0,0,1,0,1,0,1,1,0;  %8
    1,0,0,0,0,1,0,0,1]; %9

% eliminate 3
  % 1,2,3,4,5,6,7,8,9
A1=[1,0,0,1,0,0,0,0,1;  %1
    0,1,0,1,1,1,1,0,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    1,1,0,1,0,1,0,0,0;  %4
    0,1,0,0,1,0,0,1,0;  %5
    0,1,0,1,0,1,1,0,1;  %6
    0,1,0,0,0,1,1,1,0;  %7
    0,0,0,0,1,0,1,1,0;  %8
    1,0,0,0,0,1,0,0,1]; %9

% eliminate 1, add fill-in link between 4 and 9
  % 1,2,3,4,5,6,7,8,9
A2=[0,0,0,0,0,0,0,0,0;  %1
    0,1,0,1,1,1,1,0,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    0,1,0,1,0,1,0,0,1;  %4
    0,1,0,0,1,0,0,1,0;  %5
    0,1,0,1,0,1,1,0,1;  %6
    0,1,0,0,0,1,1,1,0;  %7
    0,0,0,0,1,0,1,1,0;  %8
    0,0,0,1,0,1,0,0,1]; %9

% eliminate 5, add fill-in link between 2 and 8
  % 1,2,3,4,5,6,7,8,9
A3=[0,0,0,0,0,0,0,0,0;  %1
    0,1,0,1,0,1,1,1,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    0,1,0,1,0,1,0,0,1;  %4
    0,0,0,0,0,0,0,0,0;  %5
    0,1,0,1,0,1,1,0,1;  %6
    0,1,0,0,0,1,1,1,0;  %7
    0,1,0,0,0,0,1,1,0;  %8
    0,0,0,1,0,1,0,0,1]; %9

% eliminate 8
  % 1,2,3,4,5,6,7,8,9
A4=[0,0,0,0,0,0,0,0,0;  %1
    0,1,0,1,0,1,1,0,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    0,1,0,1,0,1,0,0,1;  %4
    0,0,0,0,0,0,0,0,0;  %5
    0,1,0,1,0,1,1,0,1;  %6
    0,1,0,0,0,1,1,0,0;  %7
    0,0,0,0,0,0,0,0,0;  %8
    0,0,0,1,0,1,0,0,1]; %9

% eliminate 7
  % 1,2,3,4,5,6,7,8,9
A5=[0,0,0,0,0,0,0,0,0;  %1
    0,1,0,1,0,1,0,0,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    0,1,0,1,0,1,0,0,1;  %4
    0,0,0,0,0,0,0,0,0;  %5
    0,1,0,1,0,1,0,0,1;  %6
    0,0,0,0,0,0,0,0,0;  %7
    0,0,0,0,0,0,0,0,0;  %8
    0,0,0,1,0,1,0,0,1]; %9

% eliminate 2
  % 1,2,3,4,5,6,7,8,9
A6=[0,0,0,0,0,0,0,0,0;  %1
    0,0,0,0,0,0,0,0,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    0,0,0,1,0,1,0,0,1;  %4
    0,0,0,0,0,0,0,0,0;  %5
    0,0,0,1,0,1,0,0,1;  %6
    0,0,0,0,0,0,0,0,0;  %7
    0,0,0,0,0,0,0,0,0;  %8
    0,0,0,1,0,1,0,0,1]; %9

% eliminate 4
  % 1,2,3,4,5,6,7,8,9
A7=[0,0,0,0,0,0,0,0,0;  %1
    0,0,0,0,0,0,0,0,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    0,0,0,0,0,0,0,0,0;  %4
    0,0,0,0,0,0,0,0,0;  %5
    0,0,0,0,0,1,0,0,1;  %6
    0,0,0,0,0,0,0,0,0;  %7
    0,0,0,0,0,0,0,0,0;  %8
    0,0,0,0,0,1,0,0,1]; %9

% eliminate 6
  % 1,2,3,4,5,6,7,8,9
A8=[0,0,0,0,0,0,0,0,0;  %1
    0,0,0,0,0,0,0,0,0;  %2
    0,0,0,0,0,0,0,0,0;  %3
    0,0,0,0,0,0,0,0,0;  %4
    0,0,0,0,0,0,0,0,0;  %5
    0,0,0,0,0,0,0,0,0;  %6
    0,0,0,0,0,0,0,0,0;  %7
    0,0,0,0,0,0,0,0,0;  %8
    0,0,0,0,0,0,0,0,1]; %9


G = graph(A,'OmitSelfLoops');
H = plot(G, 'LineWidth', 2.5, 'markers', 10, 'Layout', 'circle')
title('Undirected Graph of A')
nl = H.NodeLabel;
H.NodeLabel = '';
xd = get(H, 'XData');
yd = get(H, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')

hold on
figure

G1 = graph(A1,'OmitSelfLoops');
subplot(2,4,1)
H1 = plot(G1, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 3')
nl = H1.NodeLabel;
H1.NodeLabel = '';
xd = get(H1, 'XData');
yd = get(H1, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H1,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')

G2 = graph(A2,'OmitSelfLoops');
subplot(2,4,2)
H2 = plot(G2, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 1')
nl = H2.NodeLabel;
H2.NodeLabel = '';
xd = get(H2, 'XData');
yd = get(H2, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H2,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')
highlight(H2,[4 9], 'Edgecolor', 'm');

G3 = graph(A3,'OmitSelfLoops');
subplot(2,4,3)
H3 = plot(G3, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 5')
nl = H3.NodeLabel;
H3.NodeLabel = '';
xd = get(H3, 'XData');
yd = get(H3, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H3,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')
highlight(H3,[2 8], 'Edgecolor', 'm');

G4 = graph(A4,'OmitSelfLoops');
subplot(2,4,4)
H4 = plot(G4, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 8')
nl = H4.NodeLabel;
H4.NodeLabel = '';
xd = get(H4, 'XData');
yd = get(H4, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H4,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')

G5 = graph(A5,'OmitSelfLoops');
subplot(2,4,5)
H5 = plot(G5, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 7')
nl = H5.NodeLabel;
H5.NodeLabel = '';
xd = get(H5, 'XData');
yd = get(H5, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H5,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')

G6 = graph(A6,'OmitSelfLoops');
subplot(2,4,6)
H6 = plot(G6, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 2')
nl = H6.NodeLabel;
H6.NodeLabel = '';
xd = get(H6, 'XData');
yd = get(H6, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H6,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')

G7 = graph(A7,'OmitSelfLoops');
subplot(2,4,7)
H7 = plot(G7, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 4')
nl = H7.NodeLabel;
H7.NodeLabel = '';
xd = get(H7, 'XData');
yd = get(H7, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H7,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')

G8 = graph(A8,'OmitSelfLoops');
subplot(2,4,8)
H8 = plot(G8, 'LineWidth', 2.5, 'Layout', 'circle')
title('Eliminate Node 6')
nl = H8.NodeLabel;
H8.NodeLabel = '';
xd = get(H8, 'XData');
yd = get(H8, 'YData');
text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')
highlight(H8,[1,2,3,4,5,6,7,8,9], 'NodeColor', 'red')






