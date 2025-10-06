<!-- ProfileCard.vue -->
<template>
  <div class="card shadow-sm border-0 bg-white mb-3">
    <!-- 头像区域 -->
    <div class="text-center px-4 pt-4 pb-0">
      <div class="bg-light rounded-circle d-inline-flex align-items-center justify-content-center"
        style="width: 80px; height: 80px;">
        <span class="text-muted">头像</span>
      </div>
    </div>

    <!-- 个人信息区域 -->
    <div class="card-body p-4 text-center typography-body">
      <h3 class="card-title mb-1 fw-bold">zorrooz</h3>
      <p class="card-text text-muted mb-4">
        开发者
      </p>

      <!-- 统计数据 -->
      <div class="row g-0 text-center">
        <div class="col border-end">
          <div class="fw-bold">128</div>
          <div class="text-muted">文章</div>
        </div>
        <div class="col border-end">
          <div class="fw-bold">42</div>
          <div class="text-muted">标签</div>
        </div>
        <div class="col">
          <div class="fw-bold">999K</div>
          <div class="text-muted">字数</div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'ProfileCard'
}
</script>

<style scoped></style>
