<!-- DocumentList.vue -->
<template>
  <div>
    <!-- 卡片列表 -->
    <div class="row g-3">
      <div v-for="doc in currentPageDocs" :key="doc.id" class="col-12">
        <div class="d-flex flex-column flex-lg-row border rounded shadow-sm overflow-hidden h-100">

          <!-- 左侧元信息区 -->
          <div
            class="col-12 col-lg-4 bg-white text-center text-lg-start p-3 d-flex flex-column align-items-center align-items-lg-start">
            <!-- 发布/更新日期 -->
            <div class="mb-3 text-dark">
              <div class="text-muted small">{{ monthYear(doc.date) }}</div>
            </div>

            <!-- 三级分类（面包屑式） -->
            <div class="mb-3 text-start w-100">
              <small class="text-muted">
                <template v-for="(cat, i) in doc.category" :key="i">
                  <span>{{ cat }}</span>
                  <span v-if="i < doc.category.length - 1" class="mx-1">›</span>
                </template>
              </small>
            </div>

            <!-- 标签 -->
            <div class="d-flex flex-wrap gap-1 justify-content-center justify-content-lg-start">
              <span v-for="tag in doc.tags" :key="tag" class="badge bg-light text-dark"
                style="font-size: 0.75rem; padding: 0.25em 0.4em; font-weight: normal; border: none;">
                # {{ tag }}
              </span>
            </div>
          </div>

          <!-- 右侧内容区 -->
          <div class="col-12 col-lg-8 p-3 d-flex flex-column bg-light">
            <!-- 标题 -->
            <h5 class="fw-bold text-dark mb-2 line-clamp-1">
              {{ doc.title }}
            </h5>

            <!-- 内容预览 -->
            <p class="text-secondary mb-0 flex-grow-1 line-clamp-2">
              {{ doc.description }}
            </p>

            <!-- 查看更多 -->
            <small class="text-primary mt-2">
              <a href="#" class="text-decoration-none">查看详情 →</a>
            </small>
          </div>
        </div>
      </div>
    </div>

    <!-- 分页 -->
    <nav v-if="totalPages > 1" class="d-flex justify-content-center mt-4">
      <ul class="pagination pagination-sm">
        <li class="page-item" :class="{ disabled: currentPage === 1 }">
          <a class="page-link" href="#" @click.prevent="prevPage">&laquo;</a>
        </li>
        <li v-for="page in totalPages" :key="page" class="page-item" :class="{ active: page === currentPage }">
          <a class="page-link" href="#" @click.prevent="goToPage(page)">
            {{ page }}
          </a>
        </li>
        <li class="page-item" :class="{ disabled: currentPage === totalPages }">
          <a class="page-link" href="#" @click.prevent="nextPage">&raquo;</a>
        </li>
      </ul>
    </nav>

    <!-- 无数据 -->
    <p v-if="docs.length === 0" class="text-center text-muted mt-4">
      暂无文档
    </p>
  </div>
</template>

<script>
export default {
  name: 'DocumentCardList',
  props: {
    docs: {
      type: Array,
      required: true
    },
    perPage: {
      type: Number,
      default: 5
    }
  },
  data() {
    return {
      currentPage: 1
    };
  },
  computed: {
    totalPages() {
      return Math.ceil(this.docs.length / this.perPage);
    },
    currentPageDocs() {
      const start = (this.currentPage - 1) * this.perPage;
      return this.docs.slice(start, start + this.perPage);
    }
  },
  methods: {
    prevPage() {
      if (this.currentPage > 1) this.currentPage--;
    },
    nextPage() {
      if (this.currentPage < this.totalPages) this.currentPage++;
    },
    goToPage(page) {
      if (page >= 1 && page <= this.totalPages) {
        this.currentPage = page;
      }
    },
    // 拆分日期：日 / 月+年
    day(dateStr) {
      const date = new Date(dateStr);
      return isNaN(date.getTime()) ? '' : date.getDate();
    },
    monthYear(dateStr) {
      const date = new Date(dateStr);
      if (isNaN(date.getTime())) return '';
      return date.toLocaleDateString('zh-CN', { month: 'short', year: 'numeric' });
    }
  },
  watch: {
    docs() {
      this.currentPage = 1;
    }
  }
};
</script>

<style scoped>
.line-clamp-1 {
  display: -webkit-box;
  -webkit-line-clamp: 1;
  -webkit-box-orient: vertical;
  overflow: hidden;
  text-overflow: ellipsis;
}

.line-clamp-2 {
  display: -webkit-box;
  -webkit-line-clamp: 2;
  -webkit-box-orient: vertical;
  overflow: hidden;
  text-overflow: ellipsis;
}

.shadow-sm:hover {
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
  transition: all 0.2s;
}

/* 小屏下左侧区域样式 */
@media (max-width: 991.98px) {
  .bg-light {
    border-bottom: 1px solid #eee;
  }

  .p-3 {
    padding: 1rem;
  }
}
</style>
