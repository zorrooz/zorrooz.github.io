<!-- TagCloud.vue -->
<template>
  <div class="card shadow-sm border-0 bg-white" style="max-width: 280px;">
    <div class="card-header bg-white border-0 pt-3 px-3 pb-2">
      <h6 class="text-secondary m-0 small">标签</h6>
    </div>
    <div class="card-body p-3">
      <div class="d-flex flex-wrap gap-1">
        <span class="badge bg-light text-dark fw-normal py-2 px-2" v-for="tag in tags" :key="tag.id">
          # {{ tag.name }}
        </span>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'TagCloud',
  data() {
    return {
      tags: [
        { id: 1, name: 'Vue.js' },
        { id: 2, name: 'JavaScript' },
        // ...其他标签
      ]
    }
  },
  props: {
    tagData: {
      type: Array,
      default: () => []
    }
  },
  mounted() {
    if (this.tagData.length > 0) {
      this.tags = this.tagData;
    }
  }
}
</script>

<style scoped>
/* 无自定义CSS */
</style>
