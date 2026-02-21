classdef LinOpBorn3D <  LinOp
    % LinOpBorn3D : Born forward operator. Refer to [wolf1969three] for details
    %
    % :param G: FreeSpaceKernel3D (Green function operator, $\mathbf{G}$)
    % :param Gtil: Vainikkotilde (Green function operator for Sensor, $\mathbf{\tilde{G}}$)
    % :param M: LinOpSelectorPatch/LinOpIdentity
    % :param uin: 3D incident field
    % :param uem: field with f=zeros(this.sizein) (~Incident field) on the sensor
    % :param totalfield: boolean for computing the total or scattered field
    %
    %
    % **References**
    %
    % [1] Three-Dimensional Optical Diffraction Tomography with Lippmann-Schwinger Model.,
    %     IEEE Transactions on Computational Imaging (2020), T-a. Pham, E. Soubies, 
    %     A. Ayoub, J. Lim, D. Psaltis, and M. Unser.
    %
    % **Example** H=LinOpBorn3D(G,Gtil,M,uin,uem,totalfield)
    %
    % Please refer to the LinOp superclass for documentation
    %     Copyright (C)
    %      2020 T.-A. Pham thanh-an.pham@epfl.ch
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (SetAccess = public,GetAccess = public)
        %G;
        
        MGtild;
        uin;
        u;
        uem;
        totalfield = true;
        a = 1;
    end

    methods
        
        function this = LinOpBorn3D(Gtil,M,uin,uem,totalfield)
            this.name = 'LinOpBorn3D';
            this.isDifferentiable = true;
            this.sizein = Gtil.sizein;
            uem = gpuCpuConverter(uem);
            uin = gpuCpuConverter(uin);
            %G.mtf = gpuCpuConverter(G.mtf);
            
            if isa(Gtil,'FreeSpaceKernel3Dtilde') || isa(Gtil,'FreeSpaceKernel3DtildeMat') || isa(Gtil,'Vainikkotilde')
                this.sizeout = Gtil.sizeout;%Gtil directly outputs the measurements
                this.MGtild = Gtil;
                this.uem = uem;%M must be incident field at the sensors if total field is true, (can be empty otherwise)
            elseif isa(Gtil,'Gtild3DFFT')
                this.sizeout = M.sizeout;
                Gtil.G.mtf = gpuCpuConverter(Gtil.G.mtf);
                this.MGtild = M*Gtil;%Gtil outputs scattered field on a volume that includes the sensors
                this.uem = uem;
            else
                error('Gtild not correct');
            end
            
            %this.G = G;
            if ~isempty(uin)
                Crp = LinOpSelectorPatch(size(uin),...
                    1 + (size(uin) - this.sizein)/2,(size(uin) + this.sizein)/2);
                this.uin = Crp*uin;
            end
            this.norm = gather(this.MGtild.norm*max(abs(this.uin(:)))*abs(this.a));
            this.totalfield = totalfield;
            this.u = this.uin;
        end
    end
    
    methods (Access = protected)
        function y = apply_(this,x)       
            %this.u = this.uin + this.G*(this.uin.*x); %for display
            y = this.a*this.MGtild*(this.uin.*x);
            if this.totalfield
                y = y + this.uem;
            end
        end
        
        function out = applyAdjoint_(this,y)
            %out = conj(this.a)*conj(this.uin).*(this.G'*(this.MGtild'*y));
            out = real(conj(this.a)*conj(this.uin).*(this.MGtild'*y));
        end
        
    end
end
