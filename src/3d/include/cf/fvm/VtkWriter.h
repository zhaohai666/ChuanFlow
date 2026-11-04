#pragma once

#include <string>
#include <vector>

#include "cf/common/Field.h"
#include "cf/fvm/PolyMesh.h"

namespace cf
{

/// ASCII VTK unstructured-grid (.vtu) writer; results open directly in ParaView.
/// The added fields are referenced, not copied, so each one has to outlive the
/// write it was added for.
class VtkWriter
{
public:
    explicit VtkWriter(const polyMesh& mesh) : mesh_(mesh) {}

    void addCellScalar(const std::string& name, const Field<scalar>& values);
    void addCellVector(const std::string& name, const Field<Vector>& values);

    void write(const std::string& path) const;

private:
    const polyMesh& mesh_;
    std::vector<std::pair<std::string, const Field<scalar>*>> scalars_;
    std::vector<std::pair<std::string, const Field<Vector>*>> vectors_;
};

} // namespace cf
