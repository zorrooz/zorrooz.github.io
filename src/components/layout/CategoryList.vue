<!-- CategoryList.vue -->
<template>
  <div class="card shadow-sm border-0 bg-white" style="max-width: 280px;">
    <div class="card-header bg-white border-0 pt-3 px-3 pb-2">
      <h6 class="text-secondary m-0 small">分类</h6>
    </div>
    <ul class="list-group list-group-flush">
      <li class="list-group-item d-flex justify-content-between align-items-center px-3 py-2"
        v-for="category in categories" :key="category.id">
        <span class="text-dark">{{ category.name }}</span>
        <span class="badge bg-light text-dark fw-normal">{{ category.count }}</span>
      </li>
    </ul>
  </div>
</template>

<script>
export default {
  name: 'CategoryList',
  data() {
    return {
      categories: [
        { id: 1, name: 'Vue.js', count: 12 },
        { id: 2, name: 'JavaScript', count: 8 },
        { id: 3, name: 'CSS', count: 6 },
        // ...其他分类
      ]
    }
  },
  props: {
    categoryData: {
      type: Array,
      default: () => []
    }
  },
  mounted() {
    if (this.categoryData.length > 0) {
      this.categories = this.categoryData;
    }
  }
}
</script>

<style scoped>
/* 无自定义CSS */
</style>
