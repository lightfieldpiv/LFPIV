function out = resample_flow3d(uv, sz)
% function out = resample_flow(uv, factor, method)
%RESAMPLE_FLOW   Resample flow field
%   OUT = RESAMPLE_FLOW(IN, FACTOR[, METHOD]) resamples (resizes) the flow
%   field IN using a factor of FACTOR.  The optional argument METHOD
%   specifies the interpolation method ('bilinear' (default) or
%   'bicubic'). 
%  
%   This is a private member function of the class 'hs_optical_flow'. 
%
% Authors: Deqing Sun, Department of Computer Science, Brown University
%          Stefan Roth, Department of Computer Science, TU Darmstadt
% Contact: dqsun@cs.brown.edu, sroth@cs.tu-darmstadt.de
% $Date: 2008-10-28$
% $Revision: 0 $
%
% Copyright 2007-2008, Brown University, Providence, RI. USA
% 		     TU Darmstadt, Darmstadt, Germany 
% 
%                         All Rights Reserved
% 
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and Brown University not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
% 
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR BROWN UNIVERSITY BE LIABLE FOR
% ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.  

  
%   % Make bilinear the default method
%   if (nargin < 3)
%     method = 'bilinear';
%   end

%  ratio = sz(1) / size(uv,1);
    ratio = sz ./ size(uv(:,:,:,1));
%   u     = imresize(uv(:,:,:,1), sz, method)*ratio;
%   v     = imresize(uv(:,:,:,2), sz, method)*ratio;
  u     = resize(uv(:,:,:,1), sz)*ratio(1);
  v     = resize(uv(:,:,:,2), sz)*ratio(2);
  w     = resize(uv(:,:,:,3), sz)*ratio(3);
  out   = cat(4, u, v, w);